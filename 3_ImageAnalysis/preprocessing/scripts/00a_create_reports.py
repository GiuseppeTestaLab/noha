# This script processes microscopy images to automatically segment organoids, 
# generate masks and bounding boxes, and create PDF reports with visualizations 
# and instructions for manual annotation if needed.

import pandas as pd
import aicsimageio
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
import tifffile as tif
import os
import cv2

from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Image, Paragraph
from reportlab.lib.styles import getSampleStyleSheet
styles = getSampleStyleSheet()

import warnings
warnings.filterwarnings('ignore')

from stardist.models import StarDist2D
from csbdeep.utils import normalize

import argparse

# prints a list of available models
StarDist2D.from_pretrained()

# creates a pretrained model
model = StarDist2D.from_pretrained('2D_versatile_fluo')

from img_functions import * 

parser = argparse.ArgumentParser(description='Save masks, measures and thresholded images')
parser.add_argument('--input', type=str, help='Input image path')
parser.add_argument('--output', type=str, help='Output folder path')

args = parser.parse_args()

input_path = args.input
image_name = input_path.split('/')[-1].split('.')[0]

output_path = args.output

img = aicsimageio.AICSImage(input_path)
pixel_size = img.physical_pixel_sizes.X
pixel_size_mm = (pixel_size / 1000)**2 # find area of a single pixel in mm2

morphological_properties = ['label', 'area', 'bbox', 'centroid', 'eccentricity', 'equivalent_diameter', 'extent', 'major_axis_length', 'minor_axis_length', 'orientation', 'perimeter', 'solidity']
intensity_prop = ['label', 'intensity_mean', 'intensity_min', 'intensity_max']
manual_bbox_path = '/group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/'
#bbox_files = glob.glob(manual_bbox_path + f'*.csv')
#print(bbox_files)


# %%

for scene in img.scenes:
    
    print('Processing scene: ', scene)
    
    img.set_scene(scene)

    data = img.get_image_data('XY', C = 3)
    mask_width = data.shape[0]  
    mask_height = data.shape[1]
    
    if f'{image_name}_{scene}_manual_bbox.csv' not in os.listdir(manual_bbox_path):
        
        labels_tot, _ = model.predict_instances(normalize(data), n_tiles = model._guess_n_tiles(data))
        mask_org = labels_tot > 0
        mask_org = np.int64(mask_org)
        mask_org = mask_org.astype(np.uint8)
        regions = measure.regionprops(measure.label(mask_org))
        bounding_boxes = [region.bbox for region in regions]
        fused_bounding_boxes = fuse_overlapping_boxes(bounding_boxes)
        areas = [get_bbox_area(i) for i in fused_bounding_boxes] # we compute the area of each bounding box
        fused_bounding_boxes = [box for box, area in zip(fused_bounding_boxes, areas) if area > 1e6]
        fused_bounding_boxes = fuse_overlapping_boxes(fused_bounding_boxes)
        
        mask_total = np.zeros((mask_width, mask_height), dtype=np.uint8)
            
        vertices = {}
        
        fig, ax = plt.subplots(1, 1, figsize = (10,10))
        ax.imshow(data, cmap = 'gray', vmin = np.quantile(data, 0.1), vmax = np.quantile(data, 0.99))
        
        for i, box in enumerate(fused_bounding_boxes):
            
            mask = np.zeros((mask_width, mask_height), dtype=np.uint8)
            
            min_row, min_col, max_row, max_col = box

            vertices[i] = np.array([
                [min_col, min_row],  # Top-left corner
                [max_col, min_row],  # Top-right corner
                [max_col, max_row],  # Bottom-right corner
                [min_col, max_row]   # Bottom-left corner
            ])
                    
            # OpenCV expects points as integers
            pts = vertices[i].round().astype(np.int32)

            # Reshape for fillPoly: needs (1, N, 2)
            pts = pts.reshape((-1, 1, 2))

            # Fill polygon
            cv2.fillPoly(mask_total, [pts], color=1)
            cv2.fillPoly(mask, [pts], color=int(i) + 1)
            
            ax.contour(mask, alpha = .7, colors='red', linewidths = 1)
        
        tif.imwrite(f'{output_path}/masks_organoids/{image_name}_{scene}_organoid_mask.tif', mask_total)
        plt.savefig(f'{output_path}/check_organoids/{image_name}_scene_{scene}.png', dpi=300, bbox_inches='tight')
        
        axis0 = []
        axis1 = []
        indexes = []
        for i in vertices:
            for j in range(len(vertices[i])):
                indexes.append(i)
                axis0.append(vertices[i][j,0])
                axis1.append(vertices[i][j,1])
        
        df = pd.DataFrame({'index': indexes, 'axis-0': axis0, 'axis-1': axis1})
        df.to_csv(f'{manual_bbox_path}{image_name}_{scene}_manual_bbox.csv')
        
        del mask, df
    
    else:
        print(f'Found manual bbox file {image_name}_{scene}_manual_bbox.csv')
        df_vertices = pd.read_csv(f'{manual_bbox_path}{image_name}_{scene}_manual_bbox.csv')
        mask_total = np.zeros((mask_width, mask_height), dtype=np.uint8)
        vertices = {}
        
        fig, ax = plt.subplots(1, 1, figsize = (10,10))
        ax.imshow(data, cmap = 'gray', vmin = np.quantile(data, 0.1), vmax = np.quantile(data, 0.99))

        for i in df_vertices['index'].unique():
            
            sub = df_vertices[df_vertices['index'] == i]
            vertices[i] = list(zip(sub['axis-0'], sub['axis-1']))
            vertices[i] = np.array(vertices[i])
            
            mask = np.zeros((mask_width, mask_height), dtype=np.uint8)
            
            # OpenCV expects points as integers
            pts = vertices[i].round().astype(np.int32)

            # Reshape for fillPoly: needs (1, N, 2)
            pts = pts.reshape((-1, 1, 2))

            # Fill polygon
            cv2.fillPoly(mask, [pts], color=int(i) + 1)
            cv2.fillPoly(mask_total, [pts], color=int(i) + 1)
            
            ax.contour(mask, alpha = .7, colors='red', linewidths = 1)

        plt.savefig(f'{output_path}/check_organoids/{image_name}_scene_{scene}.png', dpi=300, bbox_inches='tight')
        
        tif.imwrite(f'{output_path}/masks_organoids/{image_name}_{scene}_organoid_mask.tif', mask_total)     
        
    
    report = SimpleDocTemplate(f'/group/testa/Project/EndPoints/TPSSU/analysis/complete_report/{image_name}_{scene}_report.pdf', pagesize=A4)
    elements = []
    

    overview_image = Image(output_path + '/check_organoids/' + image_name + '_scene_' + str(scene) + '.png', width=500, height=500)
    elements.append(overview_image)
        
    singularity_command = f"""
                        <br/>
                        If the automatic segmentation didn't work as expected, you can manually annotate the organoids with napari using the following command: <br/>
                        <br/>
                        
                        
                        singularity run -B /group/testa/ docker://alessiavalenti/imaging:TMA-0.0.1 napari
                        
                        <br/>
                        <br/> Save the txt output file as f'{image_name}_{scene}_manual_bbox.csv' at the following path /group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/
                        """

    elements.append(Paragraph(singularity_command, style=styles['Normal']))
    report.build(elements)




# %%
