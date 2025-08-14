import pandas as pd
import aicsimageio
import numpy as np
import tifffile as tif
import cv2

from skimage import filters, measure

import warnings
warnings.filterwarnings('ignore')

import argparse

def crop_image(img, bbox, offset = 300):
    min_row, min_col, max_row, max_col = bbox
    cropped = img[min_row - offset:max_row + offset, min_col - offset :max_col + offset]
    return cropped

import sys
sys.path.append('/group/testa/Users/alessia.valenti/NewImageAnalysis/EmbryoNet/')
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

name = input_path.split('/')[-1].split('.')[0]
chan = 'channel_2'

morphological_properties = ['label', 'area', 'bbox', 'centroid', 'eccentricity', 'equivalent_diameter', 'extent', 'major_axis_length', 'minor_axis_length', 'orientation', 'perimeter', 'solidity']
intensity_prop = ['label', 'intensity_mean', 'intensity_min', 'intensity_max']
manual_bbox_path = '/group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/'
#bbox_files = glob.glob(manual_bbox_path + f'*.csv')
#print(bbox_files)
cyto_mask_measures = {}
properties = ['area', 'label', 'intensity_mean', 'intensity_min', 'intensity_max']
# %%

for scene in img.scenes:
    
    print('Processing scene: ', scene)

    scanreg = scene
    
    img.set_scene(scene)
    
    df_vertices = pd.read_csv(f'{manual_bbox_path}{image_name}_{scene}_manual_bbox.csv')
    print(f'Found manual bbox file {image_name}_{scene}_manual_bbox.csv')
    # mask = np.zeros((mask_width, mask_height), dtype=np.uint8)
    vertices = {}

    data = img.get_image_data('XY', C = 2)
    mask_width = data.shape[0]  
    mask_height = data.shape[1] 

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
        cv2.fillPoly(mask, [pts], color=1)
                
        region = measure.regionprops(measure.label(mask))
        
        single_organoid_tif = crop_image(data, region[0].bbox, offset = 0)
        
        flattened_img = np.squeeze(single_organoid_tif).flatten()
        flattened_img_log = np.log1p(flattened_img)
    
        img_log = np.log1p(single_organoid_tif)
    
        thresholds = filters.threshold_multiotsu(flattened_img_log, classes=4)
    
        images_transformed = np.squeeze(img_log)
        cyto_mask = img_log > thresholds[1] 

        mask = crop_image(mask, region[0].bbox, offset = 0)

        cyto_mask = np.where(mask, cyto_mask, 0)

        rep = 'rep_' + str(i)
        
        cyto_mask_measures[f'{name}_{scanreg}_{rep}_{chan}'] = measure.regionprops_table(cyto_mask.astype('int'), intensity_image=single_organoid_tif, properties=properties)
        
        tif.imwrite(f"{output_path}/SOX2_signal_masks/{name}_{scanreg}_{rep}_{chan}_mask3_log1p.tif", cyto_mask)

    #print(f'No manual bbox file for scene {scene}.')
            

# %%
import pickle

with open(f"{output_path}/SOX2_signal_mask_measures/{name}_SOX2_signal_mask_measures.pkl", "wb") as f:
    pickle.dump(cyto_mask_measures, f)