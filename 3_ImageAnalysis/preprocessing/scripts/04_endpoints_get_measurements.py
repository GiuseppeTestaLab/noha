# %%
import pandas as pd

from skimage import measure

import glob
import argparse

import cv2
import tifffile as tif

#from sklearn.preprocessing import StandardScaler
#scaler = StandardScaler()

import warnings
warnings.filterwarnings('ignore')

# %%
parser = argparse.ArgumentParser(description='Save masks, measures and thresholded images')
parser.add_argument('--input', type=str, help='Input image path')
parser.add_argument('--output', type=str, help='Output folder path')

args = parser.parse_args()

input_path = args.input
numb_image = input_path.split('/')[-1].split('.')[0].split('_')[3]

output_path = args.output

# %%

#input_path = '/group/testa/Project/EndPoints/TPSSU/Pictures/20240424_manuel_lessi_0148.czi'
#numb_image = input_path.split('/')[-1].split('.')[0].split('_')[3]

manual_bbox_path = '/group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/'
#output_path = '/group/testa/Project/EndPoints/TPSSU/analysis/'

files = glob.glob('/group/testa/Project/EndPoints/TPSSU/analysis/single_tif_organoids/*.tif')
files = [i for i in files if numb_image in i]

czi_name = ['_'.join(i.split('/')[-1].split('_')[0:4]) for i in files]
scanregions = ['_'.join(i.split('/')[-1].split('_')[4:5]) for i in files]
org_rep = ['_'.join(i.split('/')[-1].split('_')[5:7]) for i in files]
channel = ['_'.join(i.split('/')[-1].split('_')[7:9]) for i in files]

morphological_properties = ['label', 'area', 'bbox', 'centroid', 'eccentricity', 'equivalent_diameter', 'extent', 'major_axis_length', 'minor_axis_length', 'orientation', 'perimeter', 'solidity']
intensity_prop = ['label', 'intensity_mean', 'intensity_min', 'intensity_max']

# %%
for file, name, scanreg, rep, chan in zip(files, czi_name, scanregions, org_rep, channel):
        
    image = tif.imread(file)
    nuclei_mask =  tif.imread(f'{output_path}/mask_nuclei/{name}_{scanreg}_{rep}_organoid_mask.tif')
    
    measures = measure.regionprops_table(nuclei_mask, image, properties = morphological_properties + intensity_prop)
    measures = pd.DataFrame(measures).set_index('label')
    
    measures.to_csv(f'{output_path}/measurements_normalized/{name}_{scanreg}_{rep}_{chan}_measures_normalized.csv')
# %%
