import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.utils import ImageReader
from io import BytesIO

import glob

import tifffile as tif

import sys
sys.path.append('/group/testa/Users/alessia.valenti/NewImageAnalysis/')
from utilities import * 

import warnings
warnings.filterwarnings('ignore')

import argparse

parser = argparse.ArgumentParser(description='Save masks, measures and thresholded images')
parser.add_argument('--input', type=str, help='Input image path')
parser.add_argument('--output', type=str, help='Output folder path')

args = parser.parse_args()

input_path = args.input
numb_image = input_path.split('/')[-1].split('.')[0].split('_')[3]

output_path = args.output

manual_bbox_path = '/group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/'

files = glob.glob('/group/testa/Project/EndPoints/TPSSU/analysis/single_tif_organoids/*.tif')
files = [i for i in files if numb_image in i]

czi_name = ['_'.join(i.split('/')[-1].split('_')[0:4]) for i in files]
scanregions = ['_'.join(i.split('/')[-1].split('_')[4:5]) for i in files]
org_rep = ['_'.join(i.split('/')[-1].split('_')[5:7]) for i in files]
channel = ['_'.join(i.split('/')[-1].split('_')[7:9]) for i in files]

morphological_properties = ['label', 'area', 'bbox', 'centroid', 'eccentricity', 'equivalent_diameter', 'extent', 'major_axis_length', 'minor_axis_length', 'orientation', 'perimeter', 'solidity']
intensity_prop = ['label', 'intensity_mean', 'intensity_min', 'intensity_max']

for file, name, scanreg, rep, chan in zip(files, czi_name, scanregions, org_rep, channel):
        
    image = tif.imread(file)
    DAPI = tif.imread(f"/group/testa/Project/EndPoints/TPSSU/analysis/single_tif_organoids/{name}_{scanreg}_{rep}_channel_1_organoid_mask.tif")
    
    pdf_buffer = BytesIO()
    c = canvas.Canvas(pdf_buffer, pagesize=letter)
    
    nuclei_mask =  tif.imread(f'{output_path}/mask_nuclei/{name}_{scanreg}_{rep}_organoid_mask.tif')
    measures_df = pd.read_csv(f'{output_path}/measurements/{name}_{scanreg}_{rep}_{chan}_measures.csv', index_col=0)

    fig = plot_nuclei_contour(image, nuclei_mask, title = 'Nuclei segmentation before filtering', crop_center=False, cmap='Greys_r', save = False, path = None)
    
    img_buffer = BytesIO()
    fig.savefig(img_buffer, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer.seek(0)

    image_c = ImageReader(img_buffer)
    c.drawImage(image_c, x=100, y=500, width=500, height=250)
    
    fig = plot_nuclei_contour(image, nuclei_mask, title = 'Nuclei segmentation before filtering - cropped', crop_center=True, cmap='Greys_r', save = False, path = None)
    
    img_buffer = BytesIO()
    fig.savefig(img_buffer, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer.seek(0)

    image_c = ImageReader(img_buffer)
    c.drawImage(image_c, x=100, y=100, width=500, height=250)
    
    c.showPage()
    
    nuclei_mask = remove_small_instances(nuclei_mask,measures_df, measures_df['area'].mean() - 1 * measures_df['area'].std())
    nuclei_mask = remove_large_instances(nuclei_mask,measures_df, measures_df['area'].mean() + 4 * measures_df['area'].std())
    
    fig = plot_nuclei_contour(image, nuclei_mask, title = 'Nuclei segmentation after filtering', crop_center=False, cmap='Greys_r', save = False, path = None)
    
    img_buffer_2 = BytesIO()
    fig.savefig(img_buffer_2, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer_2.seek(0)

    image_c_2 = ImageReader(img_buffer_2)
    c.drawImage(image_c_2, x=100, y=500, width=500, height=250)
    
    fig = plot_nuclei_contour(image, nuclei_mask, title = 'Nuclei segmentation after filtering - cropped', crop_center=True, cmap='Greys_r', save = False, path = None)
    
    img_buffer_2 = BytesIO()
    fig.savefig(img_buffer_2, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer_2.seek(0)

    image_c_2 = ImageReader(img_buffer_2)
    c.drawImage(image_c_2, x=100, y=100, width=500, height=250)
    
    c.showPage()
    
    fig = plot_otsu_histogram(measures_df['intensity_mean'], title = "Intensity histogram and Otsu's thresolds", save = False, path = None)
    
    img_buffer_2 = BytesIO()
    fig.savefig(img_buffer_2, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer_2.seek(0)

    image_c_2 = ImageReader(img_buffer_2)
    c.drawImage(image_c_2, x=100, y=500, width=500, height=250)
    
    c.showPage()
    
    fig = visualize_nuclei_across_thresholds(image, nuclei_mask, measures_df['intensity_mean'], cmap='Greys_r', crop_center=False, save = False, path = None)
    
    img_buffer_2 = BytesIO()
    fig.savefig(img_buffer_2, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer_2.seek(0)

    image_c_2 = ImageReader(img_buffer_2)
    c.drawImage(image_c_2, x=100, y=500, width=500, height=200)
    
    fig = visualize_nuclei_across_thresholds(image, nuclei_mask, measures_df['intensity_mean'], cmap='Greys_r', crop_center=True, save = False, path = None)
    
    img_buffer_2 = BytesIO()
    fig.savefig(img_buffer_2, format='PNG', bbox_inches='tight') 
    plt.close(fig)
    img_buffer_2.seek(0)

    image_c_2 = ImageReader(img_buffer_2)
    c.drawImage(image_c_2, x=100, y=100, width=500, height=200)
    
    c.save()

    with open(f"{output_path}threhsolds_report/{name}_{scanreg}_{rep}_{chan}.pdf", "wb") as f:
        f.write(pdf_buffer.getvalue())
    