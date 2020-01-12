# -*- coding: utf-8 -*-
"""
Created on Sun May 14 13:54:22 2017

@author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg

@descript:  A pipeline that performs a range of intensity quantifications
            on image stacks where the overall tissue has previously been
            masked with `quant.preprocessing`.
"""


#------------------------------------------------------------------------------

### Imports

# Common external imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import sys, os

# Other external imports
from tifffile import imread, imsave


#------------------------------------------------------------------------------

### Main function

def intensity_analysis_pipeline(fpath, fname, channels, res, lumen_region,
                                save=False, verbose=False, show=False):
    """Pipeline that performs various intensity measurements.
    
    Parameters
    ----------
    fpath : str
        Path to the data directory.
    fname : str
        Name of the input file.
    channels : list of strings
        Name of each channel in the input file.
        Will be used in output dictionary.
    res : list of floats
        Pixel size in each dimension: [z, y, x].
    lumen_region : float
        Radius of the sphere surrounding the lumen 
        that is considered the 'lumen region'.
        
    Returns
    -------
    int_dict :  dict
        Dictionary with all computed intensities.
    """
    
    #--------------------------------------------------------------------------

    ### Load data    

    # Import input image
    stack = imread(os.path.join(fpath, fname))
    stack = np.rollaxis(stack, 1)  # Roll channels to the front

    # Report
    if verbose:
        print "  Imported image", fname
        print "  Image shape is:", stack.shape

    # Slice out membranes
    lyn = stack[0, ...]

    # Show
    if show:
        plt.imshow(lyn[lyn.shape[0]//2,:,:], cmap='gray', interpolation='none')
        plt.title('raw')
        plt.axis('off')
        plt.show()
        
    # Reconstruct the overall neuromast mask
    msk = lyn != 0
    msk = ndi.binary_fill_holes(msk)
    
    # Show
    if show:
        plt.imshow(msk[msk.shape[0]//2,:,:], cmap='gray', interpolation='none')
        plt.title('mask')
        plt.axis('off')
        plt.show()

    
    #--------------------------------------------------------------------------

    ### Mask the membrane region

    # Gaussian smoothing
    lyn_smooth = ndi.gaussian_filter(lyn, sigma=3)  # Param

    # Show
    if show:
        plt.imshow(lyn_smooth[lyn_smooth.shape[0]//2,:,:], 
                   cmap='gray', interpolation='none')
        plt.title('smoothed')
        plt.axis('off')
        plt.show()
    
    # Adaptive thresholding
    lyn_bg = ndi.uniform_filter(lyn_smooth, size=(33,100,100))  # Param
    lyn_thresh = lyn_smooth > lyn_bg

    # Show
    if show:
        plt.imshow(lyn_thresh[lyn_thresh.shape[0]//2,:,:], 
                   cmap='gray', interpolation='none')
        plt.title('bgsubbed')
        plt.axis('off')
        plt.show()
          
    # Closing to smoothen outlines    
    mem = ndi.binary_closing(lyn_thresh, iterations=5)   # Param
    
    # Erosion to thin a little
    mem = ndi.binary_erosion(mem, iterations=4)   # Param
                             
    # Show
    if show:
        plt.imshow(mem[mem.shape[0]//2,:,:], 
                   cmap='gray', interpolation='none')
        plt.title('cleaned & eroded')
        plt.axis('off')
        plt.show()
    
    # Save the membrane mask
    if save:
        imsave(os.path.join(fpath,fname[:-10]+'memmask.tif'), 
               mem.astype(np.uint8), bigtiff=True)
    
    
    #--------------------------------------------------------------------------
    
    ### Mask the lumen region
    
    # Import lumen data
    lumen = 'none'
    with open(os.path.join(fpath, r"metadata.txt"),"r") as infile:
        for line in infile.readlines():
            line = line.strip()
            line = line.split('\t')
            if line[0] in fname:
                lumen = np.array([int(value) for value in line[1:4]])
                break
    if lumen is 'none':
        raise Exception("Appropriate lumen metadata not found. Aborting!")
    
    # Create 3D mask: input params
    res    = res
    radius = lumen_region

    # Create 3D mask: adjust for resolution
    gridScale    = radius * 2
    gridSize     = [int(gridScale/r)+1 if int(gridScale/r)%2!=0
                    else int(gridScale/r) 
                    for r in res]
    centerCoords = [int(gridScale/2) for gs in gridSize]
    
    # Create 3D mask: create mask
    zz,yy,xx = np.mgrid[0:gridSize[0], 0:gridSize[1], 0:gridSize[2]] 
    sum_coords = (zz*res[0] - centerCoords[0]) ** 2 + \
                 (yy*res[1] - centerCoords[1]) ** 2 + \
                 (xx*res[2] - centerCoords[2]) ** 2
    lm3D = sum_coords <= radius ** 2

    # Place 3D mask at lumen position: prep
    lumen_mask = np.zeros_like(lyn, dtype=np.bool)
    lumen_bbox = [[lumen[0]-lm3D.shape[0]//2, lumen[0]+lm3D.shape[0]//2],
                  [lumen[1]-lm3D.shape[1]//2, lumen[1]+lm3D.shape[1]//2],
                  [lumen[2]-lm3D.shape[2]//2, lumen[2]+lm3D.shape[2]//2]]
    
    # Place 3D mask at lumen position: handle awkward boundary conditions
    if lumen_bbox[0][0] < 0:
        lumen_bbox[0][0] = 0
        lm3D = lm3D[np.abs(lumen[0]-lm3D.shape[0]//2):,:,:]
    if lumen_bbox[1][0] < 0:
        lumen_bbox[1][0] = 0
        lm3D = lm3D[:,np.abs(lumen[1]-lm3D.shape[1]//2):,:]
    if lumen_bbox[2][0] < 0:
        lumen_bbox[2][0] = 0
        lm3D = lm3D[:,:,np.abs(lumen[2]-lm3D.shape[2]//2):]
    if lumen_bbox[0][1] > lumen_mask.shape[0]:
        lumen_bbox[0][1] = lumen_mask.shape[0]
        lm3D = lm3D[:lumen[0] + lm3D.shape[0]//2 - lumen_mask.shape[0],:,:]
    if lumen_bbox[1][1] > lumen_mask.shape[1]:
        lumen_bbox[1][1] = lumen_mask.shape[1]
        lm3D = lm3D[:,:lumen[1] + lm3D.shape[1]//2 - lumen_mask.shape[1],:]
    if lumen_bbox[2][1] > lumen_mask.shape[2]:
        lumen_bbox[2][1] = lumen_mask.shape[2]
        lm3D = lm3D[:,:,:lumen[2] + lm3D.shape[2]//2 - lumen_mask.shape[2]]
    
    # Place 3D mask at lumen position: paste mask
    lumen_mask[lumen_bbox[0][0] : lumen_bbox[0][1],
               lumen_bbox[1][0] : lumen_bbox[1][1],
               lumen_bbox[2][0] : lumen_bbox[2][1]] = lm3D

    # Clip with outside mask
    lumen_mask = np.logical_and(lumen_mask, msk)

    
    #--------------------------------------------------------------------------
    
    ### Mask the basal region (very simple approach...)
    
    # Erode mask a bit for clean-up
    lower_msk = ndi.binary_erosion(msk, iterations=10)

    # Find basal quarter of extent in z
    z_extent = np.where(lower_msk==1)[0]
    z_min, z_max = (z_extent.min(), z_extent.max())
    midpoint = int(z_max - (z_max-z_min)/4.0)
    
    # Create masked basal part of stack 
    lower_msk = np.zeros_like(lower_msk)
    lower_msk[midpoint:, : , :] = 1
    
    # Clip with overall neuromast mask
    lower_msk = np.logical_and(lower_msk, msk)
    
    
    #--------------------------------------------------------------------------
    
    ### Quantify intensities
    
    # Initialize output dict
    int_dict = {}  
    
    # Measurements (other channels)
    for i in range(stack.shape[0]):
        c = channels[i]
        
        int_dict["all_"+c] = np.mean(stack[i,...][msk])
        int_dict["mem_"+c] = np.mean(stack[i,...][np.logical_and(mem, ~lumen_mask)])
        int_dict["cyt_"+c] = np.mean(stack[i,...][np.logical_and(msk, ~np.logical_or(mem,lumen_mask))])  
        int_dict["lum_"+c] = np.mean(stack[i,...][lumen_mask])
        int_dict["bas_all_"+c] = np.mean(stack[i,...][lower_msk])
        int_dict["bas_mem_"+c] = np.mean(stack[i,...][np.logical_and(mem, lower_msk)])
        int_dict["bas_cyt_"+c] = np.mean(stack[i,...][np.logical_and(~mem, lower_msk)])


    #--------------------------------------------------------------------------

    ### Save results
    
    if save:
        import pickle
        outfile_path = os.path.join(fpath, fname[:-10] + 'measurements.pkl')
        with open(outfile_path, 'wb') as outfile:
            pickle.dump(int_dict, outfile, pickle.HIGHEST_PROTOCOL)
    
    
    #--------------------------------------------------------------------------

    ### Return results
    
    return int_dict


#------------------------------------------------------------------------------

### Wrapper for Multi-Processing

# Unpack parameters and forward to the main function
def intensity_analysis_pipeline_UNPACK(params):
    dirpath, fname, channels, res, lumen_region, save, verbose, show = params
    out = intensity_analysis_pipeline(dirpath, fname, channels, res, 
                                      lumen_region, save, verbose, show)
    return out


# Get filenames and run the multiprocessing
def run_multiprocessed(dirpath, fname_end, channels, res, 
                       lumen_region, processes, 
                       save=False, verbose=False, show=False):
    
    # Get filenames
    fnames = [fname for fname in os.listdir(dirpath) 
              if fname.endswith(fname_end)]
    
    # Construct parameter vector
    params = [(dirpath, fname, channels, res, lumen_region,
               save, verbose, show) for fname in fnames]   
    
    # Run multiprocessed
    from multiprocessing import pool
    currentPool = pool.Pool(processes=processes)
    output = currentPool.map(intensity_analysis_pipeline_UNPACK, params)

    # Return result
    return output


#------------------------------------------------------------------------------



