def crop_image(img, bbox, offset = 300):
    min_row, min_col, max_row, max_col = bbox
    cropped = img[min_row - offset:max_row + offset, min_col - offset :max_col + offset]
    return cropped

def crop_center_image(img, half_size = 300):

    center_y, center_x = img.shape[0] // 2, img.shape[1] // 2
    cropped_original = img[center_y-half_size:center_y+half_size, center_x-half_size:center_x+half_size]
    
    return cropped_original

def crop_random_box(image, box_width, box_height, seed = 42):
    
    from PIL import Image
    import random
    
    image = Image.fromarray(image)
    img_width, img_height = image.size
    
    random.seed(seed)

    # Ensure the box dimensions are within the image dimensions
    if box_width > img_width or box_height > img_height:
        raise ValueError("Box dimensions are larger than the image dimensions")

    # Generate random coordinates for the top-left corner of the box
    left = random.randint(0, img_width - box_width)
    top = random.randint(0, img_height - box_height)
    right = left + box_width
    bottom = top + box_height

    cropped_image = image.crop((left, top, right, bottom))
    bbox = (top, left, bottom, right)

    return cropped_image, bbox


def show_image(imgs, colors = None, n_channels = 4, physical_px_size = None, qvmax = 0.99, qvmin = 0.1):

    from microfilm.microplot import microshow
    import matplotlib.pyplot as plt
    from matplotlib_scalebar.scalebar import ScaleBar
    import copy
    import seaborn as sns

    if colors == None:

        colors = { 0: sns.dark_palette('#0072b2', as_cmap=True),
            1 : sns.dark_palette('#009e73', as_cmap=True),
            2 : sns.dark_palette('#d55e00', as_cmap=True),
            3 : sns.dark_palette('#ffffff', as_cmap = True)
        }

    colors = {i: j for i, j in zip(imgs.keys(), colors.values())}
        
    scalebars = {}
    scalebars[0] = ScaleBar(physical_px_size, units='um', color = 'white', box_alpha=0)
    for i in range(n_channels):
        scalebars[i] = copy.deepcopy(scalebars[0])

    ncols = 2
    nrows = int(n_channels / 2) 

    if ncols % 2 > 0:
        nrows += 1

    fig, ax = plt.subplots(nrows, ncols, figsize = (14 *ncols, 14 * nrows), gridspec_kw={'wspace': 0, 'hspace': 0 })
    
    fig.patch.set_facecolor('none')
    ax = ax.flatten().T

    for i, ax in zip(imgs, ax):

        microshow(
            images=imgs[i],
            cmaps=colors[i],
            unit='um', scalebar_size_in_units=150, scalebar_unit_per_pix=physical_px_size, 
            scalebar_font_size=20,
            label_text='A', label_font_size=0.04, dpi = 70, ax = ax)

        #ax.add_artist(scalebars[i])
        #ax.tick_params(left = False, right = False , labelleft = False, labelbottom = False, bottom = False)

    fig.show()
    
def remove_small_instances(instance_mask, measurements, size_threshold):
    
    import numpy as np
        
    original_shape = instance_mask.shape
    labels_areas = np.unique(instance_mask, return_counts=True)[1]

    area_label_remove = measurements.area < size_threshold
    labels_to_remove = area_label_remove[area_label_remove].index
    mask = np.in1d( instance_mask, labels_to_remove ).reshape( original_shape )
    instance_mask[mask] = 0
    
    return instance_mask

def remove_large_instances(instance_mask, measurements, size_threshold):
    
    import numpy as np
        
    original_shape = instance_mask.shape
    labels_areas = np.unique(instance_mask, return_counts=True)[1]

    area_label_remove = measurements.area > size_threshold
    labels_to_remove = area_label_remove[area_label_remove].index
    mask = np.in1d( instance_mask, labels_to_remove ).reshape( original_shape )
    instance_mask[mask] = 0
    
    return instance_mask

# Here I'm defining two functions: one to save and one to load `.pickle` files, that are files that contain Python dictionaries.

def load_object(filename):

    import pickle

    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported): ", ex)

def save_object(obj, filename):

    import pickle

    try:
        with open(filename + ".pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
            
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)



#########################à

# Here I'm preparing some utilities for the bounding box detection and fusion

def overlap(box1, box2):
    """ Check if two bounding boxes overlap. """
    min_row1, min_col1, max_row1, max_col1 = box1
    min_row1, min_col1 = min_row1 - 100, min_col1 - 100
    max_row1, max_col1 = max_row1 + 100, max_col1 + 100
    
    min_row2, min_col2, max_row2, max_col2 = box2
    min_row2, min_col2 = min_row2 - 100, min_col2 - 100
    max_row2, max_col2 = max_row2 + 100, max_col2 + 100

    # Check if there is no overlap
    if max_row1 < min_row2 or max_row2 < min_row1 or max_col1 < min_col2 or max_col2 < min_col1:
        return False
    return True

def merge_boxes(box1, box2):
    """ Merge two bounding boxes into a single bounding box. """
    min_row1, min_col1, max_row1, max_col1 = box1
    min_row2, min_col2, max_row2, max_col2 = box2

    min_row = min(min_row1, min_row2)
    min_col = min(min_col1, min_col2) 
    max_row = max(max_row1, max_row2)
    max_col = max(max_col1, max_col2)

    return (min_row, min_col, max_row, max_col)

def fuse_overlapping_boxes(bounding_boxes):
    """ Fuse all overlapping bounding boxes into a single list of non-overlapping boxes. """
    fused_boxes = []

    while bounding_boxes:
        current_box = bounding_boxes.pop(0)
        has_merged = False

        for i, other_box in enumerate(fused_boxes):
            if overlap(current_box, other_box):
                fused_boxes[i] = merge_boxes(current_box, other_box)
                has_merged = True
                break

        if not has_merged:
            fused_boxes.append(current_box)

    return fused_boxes

def get_bbox_area(bbox):
    min_row, min_col, max_row, max_col = bbox
    area = (max_row - min_row) * (max_col - min_col)
    return area

def plot_bounding_boxes(image, bounding_boxes, save = False, path = None):
    
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    
    fig, ax = plt.subplots(1, figsize=(12, 12))
    ax.imshow(image, cmap='gray')
    
    for bbox in bounding_boxes:
        min_row, min_col, max_row, max_col = bbox
        rect = Rectangle((min_col, min_row), max_col - min_col, max_row - min_row,
                         linewidth=3, edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        
    if save:
        plt.savefig(path)
    else:
        plt.show()
        
#############################################

def plot_otsu_histogram(measurements, title='Multi-Otsu threhsolding', save = False, path = None):

    """
    Plots a histogram of the given measurements with Otsu's multi-thresholds and highlights 
    regions based on the thresholds.
    Parameters:
    -----------
    measurements : array-like
        The data values for which the histogram will be plotted.
    title : str, optional
        The title of the plot. Default is 'Multi-Otsu thresholding'.
    save : bool, optional
        If True, the plot will be saved to the specified path. If False, the plot will be displayed.
        Default is False.
    path : str or None, optional
        The file path where the plot will be saved if `save` is True. If None, no file will be saved.
        Default is None.
    Returns:
    --------
    None
        This function does not return any value. It either displays or saves the plot.
    Notes:
    ------
    - The function uses Otsu's multi-thresholding to divide the data into four classes.
    - The thresholds are visualized as vertical lines on the histogram.
    - Regions beyond the thresholds are highlighted with different colors for better visualization.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns
    from skimage import filters
    import numpy as np

    # Set the style of seaborn
    sns.set(style="whitegrid")

    fig, ax = plt.subplots(1, 1, figsize = (8, 4) )    
    plt.suptitle(title)
    
    try:
        
        otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes = 4)
        
        sns.histplot(measurements, ax = ax)
        ax.axvline(otsu_thrs[0])
        ax.axvline(otsu_thrs[1])
        ax.axvline(otsu_thrs[2])
        
        ax.axvspan(otsu_thrs[2], measurements.max(), facecolor='red', alpha=0.5)
        ax.axvspan(otsu_thrs[1], measurements.max(), facecolor='g', alpha=0.5)
        ax.axvspan(otsu_thrs[0], measurements.max(), facecolor='y', alpha=0.5)
        
        plt.tight_layout()
    except ValueError:
        
        try:
            print("Could not divide the data in 4 classes. Using 3 classes instead.")
            
            otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes = 3)
        
            sns.histplot(measurements, ax = ax)
            ax.axvline(otsu_thrs[0])
            ax.axvline(otsu_thrs[1])
            
            ax.axvspan(otsu_thrs[1], measurements.max(), facecolor='g', alpha=0.5)
            ax.axvspan(otsu_thrs[0], measurements.max(), facecolor='y', alpha=0.5)
            
        except ValueError:
            try:
                print("Could not divide the data in 3 classes. Using 2 classes instead.")
                
                otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes = 2)
            
                sns.histplot(measurements, ax = ax)
                ax.axvline(otsu_thrs[0])
                
                ax.axvspan(otsu_thrs[0], measurements.max(), facecolor='y', alpha=0.5)
                
            except ValueError:
                print("Could not divide the data in 2 classes. Using 1 class instead.")
                
                otsu_thrs = [0]
                
                sns.histplot(measurements, ax = ax)        
    
    if save:
        plt.savefig(path)
    else:
        plt.show()
        return fig

##########################################################################################

def visualize_nuclei_across_thresholds(image, nuclei_mask, measurements, cmap='Greys_r', crop_center=False, save = False, path = None):
    """
    Visualize nuclei across different intensity thresholds using multi-Otsu thresholding.
    This function generates a visualization of nuclei in an image segmented into three confidence levels
    (high, medium, and low) based on intensity measurements. It uses multi-Otsu thresholding to determine
    the thresholds and overlays the segmented regions on the original image.
    Parameters:
    -----------
    image : numpy.ndarray
        The input image to visualize, typically a grayscale or single-channel image.
    nuclei_mask : numpy.ndarray
        The mask of nuclei regions corresponding to the input image.
    measurements : pandas.Series or numpy.ndarray
        Intensity measurements for each nucleus, used to determine confidence levels. Make sure that the index corresponds to the labels in the nuclei mask.
    cmap : list, optional
        A lcolormap to use for displaying the image. If not provided, default colormaps will be used.
    Returns:
    --------
    None
        Displays a matplotlib figure with three subplots:
        - High confidence nuclei overlayed on the image.
        - Medium confidence nuclei overlayed on the image.
        - Low confidence nuclei overlayed on the image.
    Notes:
    ------
    - The function uses multi-Otsu thresholding to divide the intensity measurements into four classes.
    - The highest three classes are used to define high, medium, and low confidence levels.
    """
    
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import filters
    
    try:
        otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes=4)
        
        df_intensity = pd.DataFrame(measurements, index=measurements.index, columns=['intensity_mean'])

        indexes_high_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-1]).flatten()].index.tolist()
        indexes_high_conf_2 = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-2]).flatten()].index.tolist()
        indexes_low_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-3]).flatten()].index.tolist()

        mask_high_conf = np.where(np.isin(nuclei_mask.flatten(), indexes_high_conf), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)
        mask_high_conf_2 = np.where(np.isin(nuclei_mask.flatten(), indexes_high_conf_2), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)
        mask_low_conf = np.where(np.isin(nuclei_mask.flatten(), indexes_low_conf), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)

        if crop_center:
            image = crop_center_image(image)
            mask_high_conf = crop_center_image(mask_high_conf)
            mask_high_conf_2 = crop_center_image(mask_high_conf_2)
            mask_low_conf = crop_center_image(mask_low_conf)

        fig, ax = plt.subplots(1, 3, figsize=(30, 15))

        ax[0].imshow(image,
                    cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
        ax[0].contour(mask_high_conf,
                    cmap='Reds', alpha=.5, linewidths=.5)

        ax[1].imshow(image,
                    cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
        ax[1].contour(mask_high_conf_2,
                    cmap='Reds', alpha=.5, linewidths=.5)

        ax[2].imshow(image,
                    cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
        ax[2].contour(mask_low_conf,
                    cmap='Reds', alpha=.5, linewidths=.5)

        ax[0].xaxis.set_visible(False)
        ax[0].yaxis.set_visible(False)
        ax[0].set_title('High threshold')

        ax[1].xaxis.set_visible(False)
        ax[1].yaxis.set_visible(False)
        ax[1].set_title('Medium threshold')

        ax[2].xaxis.set_visible(False)
        ax[2].yaxis.set_visible(False)
        ax[2].set_title('Low threshold')
        
    except ValueError:
        try:
            print("Could not divide the data in 4 classes. Using 3 classes instead.")
            otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes=3)
            
            df_intensity = pd.DataFrame(measurements, index=measurements.index, columns=['intensity_mean'])

            indexes_high_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-1]).flatten()].index.tolist()
            indexes_low_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-2]).flatten()].index.tolist()

            mask_high_conf = np.where(np.isin(nuclei_mask.flatten(), indexes_high_conf), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)
            mask_low_conf = np.where(np.isin(nuclei_mask.flatten(), indexes_low_conf), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)

            if crop_center:
                image = crop_center_image(image)
                mask_high_conf = crop_center_image(mask_high_conf)
                mask_low_conf = crop_center_image(mask_low_conf)

            fig, ax = plt.subplots(1, 2, figsize=(20, 10))

            ax[0].imshow(image,
                        cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
            ax[0].contour(mask_high_conf,
                        colors='red', alpha=.5, linewidths=.5)

            ax[1].imshow(image,
                        cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
            ax[1].contour(mask_low_conf,
                        colors='red', alpha=.5, linewidths=.5)

            ax[0].xaxis.set_visible(False)
            ax[0].yaxis.set_visible(False)
            ax[0].set_title('High threshold')

            ax[1].xaxis.set_visible(False)
            ax[1].yaxis.set_visible(False)
            ax[1].set_title('Low threshold')
            
        except ValueError:
            
            try:
                print("Could not divide the data in 3 classes. Using 2 classes instead.")
                otsu_thrs = filters.threshold_multiotsu(np.array(measurements), classes=2)
                
                df_intensity = pd.DataFrame(measurements, index=measurements.index, columns=['intensity_mean'])

                indexes_high_conf = df_intensity.iloc[np.argwhere(df_intensity['intensity_mean'] > otsu_thrs[-1]).flatten()].index.tolist()

                mask_high_conf = np.where(np.isin(nuclei_mask.flatten(), indexes_high_conf), nuclei_mask.flatten(), 0).reshape(nuclei_mask.shape)

                if crop_center:
                    image = crop_center_image(image)
                    mask_high_conf = crop_center_image(mask_high_conf)

                fig, ax = plt.subplots(1, 1, figsize=(10, 10))

                ax.imshow(image,
                            cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))
                ax.contour(mask_high_conf,
                            colors='red', alpha=.5, linewidths=.5)

                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
                ax.set_title('High threshold')

                
            except ValueError:
                print("Could not divide the data in 2 classes. Using 1 class instead.")
                otsu_thrs = [0]
                
                if crop_center:
                    image = crop_center_image(image)

                fig, ax = plt.subplots(1, 1, figsize=(10, 10))

                ax.imshow(image,
                            cmap=cmap, vmin=np.quantile(image, .05), vmax=np.quantile(image, .99))

                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
        


    if save:
        plt.savefig(path)
        plt.close(fig)
    else:
        plt.show()
        return fig
    
##########################################################################################

def plot_nuclei_contour(image, nuclei_mask, title = 'Nuclei segmentation', crop_center=False, cmap='Greys_r', save = False, path = None):

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    if crop_center:
        image = crop_center_image(image)
        nuclei_mask = crop_center_image(nuclei_mask)
        
    fig, ax = plt.subplots(1,2,figsize = (10,15))
    ax[0].imshow(
            image,
            cmap=cmap, vmin=np.quantile(image, .05), vmax = np.quantile(image, .99)
    )

    ax[1].imshow(image,
            cmap=cmap, vmin=np.quantile(image, .05), vmax = np.quantile(image, .99))

    ax[1].contour(nuclei_mask,
            colors='red', alpha = .5, linewidths = 1)

    ax[0].xaxis.set_visible(False)
    ax[0].yaxis.set_visible(False)
    ax[0].set_title(title)

    ax[1].xaxis.set_visible(False)
    ax[1].yaxis.set_visible(False)
    ax[1].set_title(title)

    if save:
        plt.savefig(path)
        plt.close(fig)
        
    else:
        plt.show()
        return fig