# -*- coding: utf-8 -*-
"""
Created on Sun May 14 13:54:22 2017

@author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg

@descript:  Preprocessing before landmark extraction for analysis of cxcr7
            distribution. This includes a masking function to mask the entire
            neuromast and a 'membrane subtraction' function that attempts to
            remove the membrane-localized signal in the cxcr7 channel based on
            the Lyn channel.
"""


#------------------------------------------------------------------------------

### Imports

# Common external imports
from __future__ import division
import os, sys
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

# Other external imports
from tifffile import imread, imsave


#------------------------------------------------------------------------------

### Prep

# Function to load data
def load_data(fname, fpath, verbose=True):
   
    # Load data
    stack = imread(os.path.join(fpath, fname))
    stack = np.rollaxis(stack,1)  # Roll channels to the front
    
    # Report
    if verbose:
        print "  Imported image", fname
        print "  Image shape is:", stack.shape

    # Return result
    return stack


#------------------------------------------------------------------------------

### Masking

# Function to create neurmast mask
def create_nm_mask(stack, show=False,
                   struct_info=([5,5,5],[2,2,2],2)):

    
    ### MEDIAN SMOOTHING
    
    # Small median filter to clean up noise
    stack = ndi.median_filter(stack, size=3)  # Param
    
    
    ### THRESHOLDING
    
    # Background percentile threshold
    # Note: The idea behind this is that we allow 5% of background pixels to
    #       be accepted as foreground, which will then be taken care of by the
    #       morphological cleaning below. For very low background values, this
    #       is a bit imprecise due to the uint format, but in general it seems
    #       to work well!
    mask = stack > np.percentile(stack[:,:50,:], 95)  # Param
    
    # Show
    if show:
        plt.imshow(mask[50,...], cmap='gray', interpolation='none')
        plt.title("Thresholded")
        plt.show()
        

    ### CLEAN-UP WITH MORPHOLOGY

    # Create 'spherical' structural element
    gridShape, centerCoords,radius = struct_info
    xx,yy,zz = np.mgrid[0:gridShape[0], 0:gridShape[1], 0:gridShape[2]] 
    sum_coords = (xx - centerCoords[0]) ** 2 + (yy - centerCoords[1]) ** 2 + (zz - centerCoords[2]) ** 2
    struct = sum_coords <= radius ** 2

    # Show struct
    if show:
        plt.imshow(struct[:,:,2], interpolation='none')
        plt.title("Structural element")
        plt.show()
    
    # Fill holes
    mask = ndi.binary_fill_holes(mask)    

    # Erode to avoid connected bg spots and to smoothen surface
    mask = ndi.binary_erosion(mask, structure=struct, iterations=6)  # Param
      
    # Retain only largest object
    labeled_mask = ndi.label(mask)[0]
    labels,counts = np.unique(labeled_mask,return_counts=True)
    counts[0] = 0
    mask = labeled_mask == labels[np.argmax(counts)]    
       
    # Dilate to reverse erosion and smoothen surface
    mask = ndi.binary_dilation(mask, structure=struct, iterations=6)  # Param

    # Show
    if show:
        plt.imshow(mask[50,...], cmap='gray', interpolation='none')
        plt.show()
        
    # Return the final mask
    return mask


#------------------------------------------------------------------------------

### Membrane Subtraction

# Function to normalize and subtract membrane intensity
def subtract_mem(stack,mem_stack):
    
    # Normalize membrane using the means
    mem_mean = np.mean(np.ma.array(mem_stack, mask=mem_stack==0))
    cx7_mean = np.mean(np.ma.array(stack, mask=stack==0))
    mem_norm = mem_stack.astype(np.float) / mem_mean * cx7_mean
                    
    # Clean up in case some values get too big
    mem_norm[mem_norm>255] = 255
    mem_norm = mem_norm.astype(np.uint8)
    
    # Subtract normed membrane
    subtracted = stack - mem_norm
    
    # Clean up negative values
    subtracted[stack < mem_norm] = 0

    # Return the result
    return subtracted


#------------------------------------------------------------------------------

### Pipeline wrapper for multiprocessing

def batch_multiprocessed(fname):
    
    # Unpack input
    fname, fpath = fname
    
    # Report
    print "## Starting run on", fname
    
    # Load data
    stack = load_data(fname, fpath, verbose=False)

    # Create neuromast mask
    mask = create_nm_mask(stack[0,...], show=False)
    
    # Apply mask
    masked_stack = np.copy(stack)
    for c in range(masked_stack.shape[0]):
        masked_stack[c,~mask] = 0 

    # Save output
    imsave(os.path.join(fpath, fname[:-4]+'_masked.tif'), 
           np.rollaxis(masked_stack,1), bigtiff=True)    
    
    # Mean-normalize and subtract membrane channel from cxcr7 channel
    memsubbed = np.zeros_like(masked_stack)
    memsubbed[0,...] = np.copy(masked_stack[0,...])
    for c in range(1, masked_stack.shape[0]):
        memsubbed[c,...] = subtract_mem(masked_stack[c,...], 
                                        masked_stack[0,...])    

    # Save output
    imsave(os.path.join(fpath, fname[:-4]+'_masked_memsub.tif'),
           np.rollaxis(memsubbed,1), bigtiff=True)


#------------------------------------------------------------------------------

### Execution wrapper for multiprocessing

def run_multiprocessed(fnames, processes):

    # Run
    from multiprocessing import pool
    currentPool = pool.Pool(processes=processes)
    currentPool.map(batch_multiprocessed, fnames)


#------------------------------------------------------------------------------



