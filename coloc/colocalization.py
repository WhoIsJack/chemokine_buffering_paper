# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:12:44 2018

@author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg

@descript:  Functions for background subtraction, region masking and intensity
            measurements for the purpose of measuring colocalization of one 
            marker within the compartments of another.
"""


#------------------------------------------------------------------------------

### Imports

# Common external imports
from __future__ import division
import os, sys
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

# Specific external imports
from skimage.filters import threshold_otsu
from scipy.stats import linregress, pearsonr


#------------------------------------------------------------------------------

### Preprocessing

# Non-adaptive global background subtraction based on background regions
def bgsub_global(stack):
    """Function for non-adaptive global background subtraction based on 
    background regions.
    
    The background is determined as the mean intensity of stripes at the
    top and bottom (in y) of the stack. The stripes are 1/10 of the image's 
    size in y and span the entire x [and z] dimension.
    
    IN:  stack; array, 1-channel image of shape ([Z],Y,X)
    OUT: bgsub; array, background subtracted image of same shape
    """
    
    # Generate background mask (stripe at top and bottom)
    bg_mask = np.zeros_like(stack, dtype=np.bool)
    bg_mask[..., :stack.shape[-2]//10 , :] = 1
    bg_mask[..., -stack.shape[-2]//10:, :] = 1

    # Compute and subtract background
    bg = np.mean(stack[bg_mask])
    bgsub = stack - bg
    bgsub[stack < bg] = 0
    
    # Done
    return bgsub


# Gaussian-adaptive local background subtraction
def bgsub_local(stack, sigma=10):
    """Function for Gaussian-adaptive background subtraction.
    
    The background is constructed by 3D Gaussian smoothing of the
    input stack with the given sigma.
    
    IN:  stack; array, 1-channel image of shape ([Z],Y,X)
         sigma; int, optional, default 10
    OUT: bgsub; array, background subtracted image of same shape
    """
    
    # Generate background by Gaussian smoothing
    bg = ndi.gaussian_filter(stack, sigma)
    
    # Background subtraction
    bgsub = stack - bg
    bgsub[stack < bg] = 0
    
    # Done
    return bgsub
    
 
#------------------------------------------------------------------------------

### Thresholding & Measurements

# Thresholding based on threshold detection
def thresh_detect(mask_stack, measure_stack):
    """Function for thresholding mask_stack based on Otsu threshold detection 
    and extracting relevant measures from measure_stack
    
    IN:  mask_stack; 8bit array, 1-channel image of shape ([Z],Y,X)
         measure_stack; array, 1-channel image of same shape as mask_stack
    OUT: threshold; int, threshold used 
         mean; float, measured mean intensity in mask
         msum; float, measured sum intensity in mask
         mean_ratio; float, mean in mask / mean in total
         msum_ratio; float, sum in mask / sum in total
    """
    
    # Get threshold
    threshold = threshold_otsu(mask_stack)
    
    # Threshold and measure
    mask = mask_stack > threshold
    mean = np.mean(measure_stack[mask])
    msum = np.sum(measure_stack[mask])
    mean_ratio = mean / np.mean(measure_stack[measure_stack>0])
    msum_ratio = msum / np.sum(measure_stack)
    
    # Return result
    return threshold, mean, msum, mean_ratio, msum_ratio #, r_otsu, r_ref, r_ratio


# Thresholding with threshold series
def thresh_series(mask_stack, measure_stack,
                  start=0, stop=256, step=5):
    """Function for performing a threshold series on mask_stack and extracting
    mean foreground intensity of measure_stack at every threshold.
    
    IN:  mask_stack; 8bit array, 1-channel image of shape ([Z],Y,X)
         measure_stack; array, 1-channel image of same shape as mask_stack
         start, stop, step; ints, optional, default 0, 256, 5
    OUT: thresholds; array, 1d-array of thresholds used
         mean_series; array, 1d-array of measured mean intensities 
         sum_series; array, 1d-array of measured sum intensities
         means_slope; linear slope of means over thresholds
         sums_slope; linear slope of sums over thresholds
    """
    
    # Get thresholds
    thresholds = np.arange(start, stop, step)
    
    # For each threshold...
    mean_series = np.zeros_like(thresholds, dtype=np.float)
    sum_series  = np.zeros_like(thresholds, dtype=np.float)
    for i,thresh in enumerate(thresholds):
        
        # Threshold and measure
        mask = mask_stack > thresh
        if np.sum(mask) > 0:
            mean = np.mean(measure_stack[mask])
            msum = np.sum(measure_stack[mask])
        else:
            mean = np.nan
            msum = np.nan
        
        # Add to mean series
        mean_series[i] = mean
        sum_series[i] = msum
        
    # Compute slopes
    means_slope = linregress(thresholds[~np.isnan(mean_series)], 
                             mean_series[~np.isnan(mean_series)])[0]
    sums_slope  = linregress(thresholds[~np.isnan(sum_series)], 
                             sum_series[~np.isnan(sum_series)])[0]
    
    # Return result
    return thresholds, mean_series, sum_series, means_slope, sums_slope


#------------------------------------------------------------------------------



