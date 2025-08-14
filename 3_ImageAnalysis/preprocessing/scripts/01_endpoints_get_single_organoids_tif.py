import pandas as pd
import aicsimageio
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif
import cv2

from skimage import measure

import warnings
warnings.filterwarnings('ignore')

import argparse

def crop_image(img, bbox, offset = 300):
    min_row, min_col, max_row, max_col = bbox
    cropped = img[min_row - offset:max_row + offset, min_col - offset :max_col + offset]
    return cropped

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


# %%

for scene in img.scenes:
    
    print('Processing scene: ', scene)
    
    img.set_scene(scene)
    
    df_vertices = pd.read_csv(f'{manual_bbox_path}{image_name}_{scene}_manual_bbox.csv')
    print(f'Found manual bbox file {image_name}_{scene}_manual_bbox.csv')
    # mask = np.zeros((mask_width, mask_height), dtype=np.uint8)
    vertices = {}
    
    for channel in range(0, 4):

        data = img.get_image_data('XY', C = channel)
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
            #print(region)
            #print(region[0].bbox)

            tif.imwrite(f'{output_path}/masks_organoids/{image_name}_{scene}_rep_{i}_organoid_mask.tif', mask) 
            
            # data = np.where(mask == 1, data, 0)
            single_organoid_tif = crop_image(data, region[0].bbox, offset = 0)
            
            tif.imwrite(f'{output_path}/single_tif_organoids/{image_name}_{scene}_rep_{i}_channel_{channel}_organoid_mask.tif', single_organoid_tif) 

            #plt.savefig(f'{output_path}/check_organoids/{image_name}_scene_{scene}.png', dpi=300, bbox_inches='tight')


    print(f'No manual bbox file for scene {scene}.')
            


# %%
