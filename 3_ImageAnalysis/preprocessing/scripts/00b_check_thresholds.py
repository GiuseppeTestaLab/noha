import pandas as pd
import aicsimageio
import matplotlib.pyplot as plt
import numpy as np
from skimage import filters, measure

import microfilm as mf
from microfilm.microplot import microshow
import os

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

manual_bbox_path = '/group/testa/Project/EndPoints/TPSSU/analysis/manual_bbox/'

if not os.path.exists(output_path + '/check_organoids/'):
    os.makedirs(output_path + '/check_organoids/')
    
if not os.path.exists(output_path + '/measurements/'):
    os.makedirs(output_path + '/measurements/')
    
if not os.path.exists(output_path + '/mask_check/'):
    os.makedirs(output_path + '/mask_check/')
    

img = aicsimageio.AICSImage(input_path)
pixel_size = img.physical_pixel_sizes.X
pixel_size_mm = (pixel_size / 1000)**2 # find area of a single pixel in mm2

morphological_properties = ['label', 'area', 'bbox', 'centroid', 'eccentricity', 'equivalent_diameter', 'extent', 'major_axis_length', 'minor_axis_length', 'orientation', 'perimeter', 'solidity']
intensity_prop = ['label', 'intensity_mean', 'intensity_min', 'intensity_max']

for scene in img.scenes:
    
    print('Processing scene: ', scene)
    
    img.set_scene(scene)
        
    condition_imgs = {}
    
    condition_imgs['channel1'] = img.get_image_data()[0][0][0]
    condition_imgs['channel2'] = img.get_image_data()[0][1][0]
    condition_imgs['channel3'] = img.get_image_data()[0][2][0]
    condition_imgs['channel4'] = img.get_image_data()[0][3][0]
    
    labels_tot, _ = model.predict_instances(normalize(condition_imgs['channel4']), n_tiles = model._guess_n_tiles(condition_imgs['channel4']))
    
    if f'{image_name}_{scene}_manual_bbox.csv' in os.listdir(manual_bbox_path):

        load_manual_annotation = pd.read_csv(f'{image_name}_{scene}_manual_bbox.csv', index_col = 0)

        fused_bounding_boxes = []
        for i in load_manual_annotation.index.unique():
            min_row, max_row = load_manual_annotation.loc[i]['axis-0'].unique()
            min_col, max_col = load_manual_annotation.loc[i]['axis-1'].unique()
            fused_bounding_boxes.append((min_row, min_col, max_row, max_col))
    
    else:
        
        mask_org = labels_tot > 0
        mask_org = np.int64(mask_org)
        mask_org = mask_org.astype(np.uint8)
        regions = measure.regionprops(measure.label(mask_org))
        bounding_boxes = [region.bbox for region in regions]
        fused_bounding_boxes = fuse_overlapping_boxes(bounding_boxes)
        areas = [get_bbox_area(i) for i in fused_bounding_boxes] # we compute the area of each bounding box
        fused_bounding_boxes = [box for box, area in zip(fused_bounding_boxes, areas) if area > 1e6]
        fused_bounding_boxes = fuse_overlapping_boxes(fused_bounding_boxes)
        
    report = SimpleDocTemplate(f'/group/testa/Project/EndPoints/TPSSU/analysis/complete_report/{image_name}_{scene}_report.pdf', pagesize=A4)
    elements = []
    
    plot_bounding_boxes(mask_org, fused_bounding_boxes, save = True, path = output_path + '/check_organoids/' + image_name + '_scene_' + str(scene) + '.png')
    
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
    
    print('Detected ', len(fused_bounding_boxes), ' organoids')
    
    measurement = {}
    new_labels = {}
    
    for i, box in enumerate(fused_bounding_boxes):
        
        measurement[f'rep{i}'] = {}
        measurement[f'rep{i}']['bbox'] = fused_bounding_boxes[i]
        measurement[f'rep{i}']['bbox_area'] = get_bbox_area(fused_bounding_boxes[i])
        
        distance_top_border = np.abs(0 - box[0])
        distance_bottom_border = np.abs(mask_org.shape[0] - box[2])

        distance_left_border = np.abs(0 - box[1])
        distance_right_border = np.abs(mask_org.shape[1] - box[3])

        all_distances = [distance_top_border, distance_bottom_border, distance_left_border, distance_right_border]
        
        if any([j < 300 for j in all_distances]) and all([j > 50 for j in all_distances]):
            offset = 50
        elif any([j < 50 for j in all_distances]):
            offset = 0
        else:
            offset = 300
        
        cropped = crop_image(condition_imgs['channel4'], fused_bounding_boxes[i], offset = offset)
        labels, _ = model.predict_instances(normalize(cropped), n_tiles = model._guess_n_tiles(cropped))
        new_labels[i] = labels
        
        measurement[f'rep{i}']['n_nuclei'] = len(np.unique(new_labels[i], return_counts = True)[0])
        
        nuclei_size = np.unique(new_labels[i], return_counts = True)[1]
        nuclei_size = nuclei_size[1:]
        
        measurement[f'rep{i}']['nuclei_areas'] = nuclei_size
        measurement[f'rep{i}']['tot_nuclei_areas'] = np.sum(measurement[f'rep{i}']['nuclei_areas'])
        measurement[f'rep{i}']['tot_nuclei_areas_mm'] = np.sum(measurement[f'rep{i}']['nuclei_areas']) * pixel_size_mm

        for channel in condition_imgs:
            
            # get cropped organoid
            cropped_channel = crop_image(condition_imgs[channel], measurement[f'rep{i}']['bbox'], offset = offset)
            measurement[f'rep{i}'][f'{channel}_img'] = cropped_channel
            
            # get intensities measures
            intensity_ch = measure.regionprops_table(label_image = new_labels[i], intensity_image = cropped_channel, properties=intensity_prop)
            measurement[f'rep{i}'][f'{channel}_measurements'] = intensity_ch
            
            # get Otsu thresholds
            otsu_thrs = filters.threshold_multiotsu(measurement[f'rep{i}'][f'{channel}_measurements']['intensity_mean'], classes = 4)
            measurement[f'rep{i}'][f'{channel}_measurements']['otsu_thresholds'] = otsu_thrs
            
            # find nuclei positive depending on threshold
            df_intensity = pd.DataFrame(measurement[f'rep{i}'][f'{channel}_measurements']['intensity_mean'], index = measurement[f'rep{i}'][f'{channel}_measurements']['label'], columns = ['intensity_mean'])
            
            indexes_high_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-1]).flatten()].index.tolist()
            indexes_high_conf_2 = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-2]).flatten()].index.tolist()
            indexes_low_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-3]).flatten()].index.tolist()

            mask_high_conf = np.where(np.isin(new_labels[i].flatten(), indexes_high_conf), new_labels[i].flatten(), 0).reshape(new_labels[i].shape)
            mask_high_conf_2 = np.where(np.isin(new_labels[i].flatten(), indexes_high_conf_2), new_labels[i].flatten(), 0).reshape(new_labels[i].shape)
            mask_low_conf = np.where(np.isin(new_labels[i].flatten(), indexes_low_conf), new_labels[i].flatten(), 0).reshape(new_labels[i].shape)

            measurement[f'rep{i}'][f'{channel}_measurements']['n_pos_thr1'] = len(np.unique(mask_high_conf))
            measurement[f'rep{i}'][f'{channel}_measurements']['n_pos_thr2'] = len(np.unique(mask_high_conf_2))
            measurement[f'rep{i}'][f'{channel}_measurements']['n_pos_thr3'] = len(np.unique(mask_low_conf))
            
            # produce thresholded images
            fig, ax = plt.subplots(3, 1, figsize = (40,15))
    
            microshow(
                        images=measurement[f'rep{i}'][f'{channel}_img'], fig_scaling=10,
                        cmaps='pure_cyan',
                        unit='um', scalebar_size_in_units=150, scalebar_unit_per_pix=pixel_size, 
                        scalebar_font_size=30,
                        label_text='A', label_font_size=0.04, dpi = 70, ax = ax[0], alpha=100)
            
            ax[0].set_title('Low confidence nuclei', fontsize = 30)
            ax[0].contour(mask_low_conf, alpha = .7, colors='red', linewidths = 0.1)
            
            microshow(
                        images=measurement[f'rep{i}'][f'{channel}_img'], fig_scaling=10,
                        cmaps='pure_cyan',
                        unit='um', scalebar_size_in_units=150, scalebar_unit_per_pix=pixel_size, 
                        scalebar_font_size=30,
                        label_text='A', label_font_size=0.04, dpi = 70, ax = ax[1], alpha=100)
            
            ax[1].set_title('Medium confidence nuclei', fontsize = 30)
            ax[1].contour(mask_high_conf_2, alpha = .7, colors='red', linewidths = 0.1)

            microshow(
                        images=measurement[f'rep{i}'][f'{channel}_img'], fig_scaling=10,
                        cmaps='pure_cyan',
                        unit='um', scalebar_size_in_units=150, scalebar_unit_per_pix=pixel_size, 
                        scalebar_font_size=30,
                        label_text='A', label_font_size=0.04, dpi = 70, ax = ax[2], alpha=100)
            
            ax[2].set_title('High confidence nuclei', fontsize = 30)
            ax[2].contour(mask_high_conf, alpha = .7, colors='red', linewidths = 0.1)
            
            # Save the thresholded images
            save_path = f'{output_path}/mask_check/{image_name}_rep{i}_{scene}_{channel}_thresholds_comparison.pdf'
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            
            #overview_image = Image(save_path, width=500, height=500)
            #elements.append(overview_image)
            
            plt.close(fig)

            morphological_properties_table = measure.regionprops_table(new_labels[i], properties=morphological_properties)
            measurement[f'rep{i}']['nuclei_morph_prop'] = morphological_properties_table

            measurement[f'rep{i}']['nuclei_morph_prop'] = pd.DataFrame(measurement[f'rep{i}']['nuclei_morph_prop'])
            
            save_object(measurement, output_path + f'/measurements/measurements_EDCs_rep{i}_{scene}_{image_name}')  
            print(f'Saved measurements for scene: {scene} at {output_path}/measurements/measurements_EDCs_{scene}_{image_name}') 
            
        