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


#########################################

# Here I'm defining a function to show all the channel images so that I don't have to copy always the same code:

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

    fig, ax = plt.subplots(nrows, ncols, figsize = (15, 7 * nrows), gridspec_kw={'wspace': 0})
    
    fig.patch.set_facecolor('none')
    ax = ax.flatten().T

    for i, ax in zip(imgs, ax):

        microshow(
            images=imgs[i], fig_scaling=5,
            cmaps=colors[i],
            unit='um', scalebar_size_in_units=150, scalebar_unit_per_pix=physical_px_size, 
            scalebar_font_size=20,
            label_text='A', label_font_size=0.04, dpi = 70, ax = ax)

        #ax.add_artist(scalebars[i])
        #ax.tick_params(left = False, right = False , labelleft = False, labelbottom = False, bottom = False)

    fig.show()

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


# Utility to crop images having a bounding box

def crop_image(img, bbox, offset = 300):
    min_row, min_col, max_row, max_col = bbox
    cropped = img[min_row - offset:max_row + offset, min_col - offset :max_col + offset]
    return cropped

