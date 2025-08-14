import tifffile as tif
import numpy as np

from skimage import filters

import glob
import sys
sys.path.append('/group/testa/Users/alessia.valenti/NewImageAnalysis/')
from utilities import * 

import argparse

parser = argparse.ArgumentParser(description='Save masks, measures and thresholded images for cytoplasmatic markers')
parser.add_argument('--input', type=str, help='Input image path')
parser.add_argument('--output', type=str, help='Output folder path')

args = parser.parse_args()


input_path = args.input
numb_image = input_path.split('/')[-1].split('.')[0].split('_')[3]

output_path = args.output

files = glob.glob('/group/testa/Project/EndPoints/TPSSU/analysis/single_tif_organoids/*.tif')
files = [i for i in files if numb_image in i and 'channel_0' in i]

czi_name = ['_'.join(i.split('/')[-1].split('_')[0:4]) for i in files]
scanregions = ['_'.join(i.split('/')[-1].split('_')[4:5]) for i in files]
org_rep = ['_'.join(i.split('/')[-1].split('_')[5:7]) for i in files]
channel = ['_'.join(i.split('/')[-1].split('_')[7:9]) for i in files]

for file, name, scanreg, rep, chan in zip(files, czi_name, scanregions, org_rep, channel):
    img = tif.imread(file)

    flattened_img = np.squeeze(img).flatten()
    flattened_img_log = np.log10(flattened_img + 1)

    img_log = np.log10(img + 1)

    thresholds = filters.threshold_multiotsu(flattened_img_log, classes=4)

    images_transformed = np.squeeze(img_log)
    mask1 = img_log > thresholds[0]
    mask2 = img_log > thresholds[1]
    mask3 = img_log > thresholds[2] 
    
    tif.imwrite(f"{output_path}/cytoplasmatic_masks/{name}_{scanreg}_{rep}_{chan}_mask3.tif", mask3)
