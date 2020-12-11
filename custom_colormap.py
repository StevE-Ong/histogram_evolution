# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

def colormap_alpha(color):
    """
    Create colormap with alpha 0-->1
    color : string 
        color: hot, Reds, viridis_r
    return : colormap with top (alpha=1) to bottom (alpha=0)
    transparent for small value
    """
    # get colormap
    ncolors = 256
    color_array = plt.get_cmap(f'{color}')(range(ncolors))
    # create alpha values
    alphas = np.zeros(ncolors)
    const = -(np.log(1e-2)-1e-2)/ncolors 
    
    for i in range(ncolors):   
        if i==0:
            alphas[i] = 0.0
        else:
            alphas[i] = 1-np.exp(-const*i+1e-2)
            
    # Define the background as white
    BG = np.asarray([1., 1., 1.,])
    
    # change alpha values
    for i in range(ncolors):
        color_array[i,:-1] = color_array[i,:-1] * alphas[i] + BG * (1.- alphas[i])   
        
    # create a colormap object
    map_object = LinearSegmentedColormap.from_list(name=f'{color}_alpha',colors=color_array)    
    # register this new colormap with matplotlib
    plt.register_cmap(cmap=map_object)

def colormap_fuse(color1,color2):
    """
    Create fused colormap of two colors
    color1 : string
        colorbar at the top
        
    color2 : string
        colorbar at the bottom
        
    return : colormap with two colors of choices 
    """
    # define top and bottom colormaps 
    top = cm.get_cmap(f'{color1}_r', 128) # r means reversed version
    bottom = cm.get_cmap(f'{color2}', 128)
    # combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
    # create a new colormaps with a name of OrangeBlue
    #new_cmp= ListedColormap(newcolors,name=f'{top}{bottom}')
    newcmp = LinearSegmentedColormap.from_list(name=f'{color1}{color2}',colors=newcolors)
    plt.register_cmap(cmap=newcmp)
    view_colormap(newcmp)

def greyscale_cmap(color):
    """
    Convert a given colormap to greyscale
    """
    cmap = cm.get_cmap(color)
    colors = cmap(np.arange(cmap.N))

    weight = [0.299, 0.587, 0.114]
    Intensity = np.sqrt(np.dot(colors[:, :3] ** 2, weight))
    colors[:, :3] =Intensity[:, np.newaxis]
        
    return LinearSegmentedColormap.from_list(cmap.name + "_grey", colors, cmap.N)

def view_colormap(color):
    """Plot a colormap with its grayscale equivalent"""
    cmap = cm.get_cmap(color)
    colors = cmap(np.arange(cmap.N))
    
    cmap = greyscale_cmap(cmap)
    greyscale = cmap(np.arange(cmap.N))
    
    fig, ax = plt.subplots(2, figsize=(4, 1),
                           subplot_kw=dict(xticks=[], yticks=[]))
    ax[0].imshow([colors], extent=[0, 10, 0, 1])
    ax[1].imshow([greyscale], extent=[0, 10, 0, 1])
# -


