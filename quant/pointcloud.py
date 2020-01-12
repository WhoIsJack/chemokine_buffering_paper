# -*- coding: utf-8 -*-
"""
Created on Sun May 14 13:54:22 2017

@author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg

@descript:  Functions for converting fluorescence intensity distributions
            into a point cloud representation and then register them to
            the image frame.
"""


#------------------------------------------------------------------------------

### Imports

# Standard external imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.ndimage as ndi
import sys, os

# Other external imports
from sklearn.decomposition import PCA


#------------------------------------------------------------------------------

### Function for landmark extraction

def generate_pcl(image, nHits, adjust_power=1, replace=False, seed=None):

    # Seed random number generator
    if seed:
        np.random.seed(seed)

    # Normalize to a total intensity of 1
    normed = image.astype(np.float) / np.sum(image)

    # Draw from distribution (without replacement)
    indx_arr = np.arange(normed.flatten().shape[0])
    hits_arr = np.random.choice(indx_arr,
                                size=nHits,
                                replace=replace,
                                p=normed.flatten())
    
    # Unravel flat index hits array
    hits = np.array(np.unravel_index(hits_arr,np.shape(normed))).T

    # Return result
    return hits

    
#------------------------------------------------------------------------------

### Main function

def point_cloud_pipeline(stack, ref_stack,
                         fpath, fname, res, num_LMs=500,
				         verbose=False, show=False):
    """Pipeline that extracts and aligns point clouds
    from intensity distributions.
    
    Parameters
    ----------
    stack : 3D numpy image
        Intensity distribution to convert to point cloud.
    ref_stack : 3D numpy image
        Intensity distribution reflecting overall tissue
        shape (usually the membranes). Used for aligning
        the cloud to the image frame and normalizing z.
    fpath : string
        Path of the source image file corresponding to the
        input stack. Used to find matching metadata.
    fname : list of strings
        File name of the source image file corresponding
        to the input stack. Used to find matching metadata.
    res : list of floats
        Pixel size in each dimension: [z, y, x].
        
    Returns
    -------
    lms : numpy array of shape (num_LMs, 3)
        Landmark coordinates in the image space (zyx).
    lms_tf : numpy array of shape (num_LMs, 3)
        Aligned and z-normalized landmark coordinates(zyx).
    lum_dist_lms : numpy array of shape (num_LMs)
        Euclidean distance of landmarks to the lumen.
    """
    
    #--------------------------------------------------------------------------

    ### Run landmark assignment

    # Run landmark assignment
    lms     = generate_pcl(stack, num_LMs, seed=42)
    ref_lms = generate_pcl(ref_stack, num_LMs, seed=42)

    # Change from pixels to um
    lms     = lms     * np.array(res)
    ref_lms = ref_lms * np.array(res)
    
    # Plot results
    if show:
        plt.scatter(lms[:,2], lms[:,1], c=lms[:,0], cmap='viridis')
        plt.title('Channel landmarks in image frame')
        plt.show()
        plt.scatter(ref_lms[:,2], ref_lms[:,1], c=ref_lms[:,0], cmap='viridis')
        plt.title('Reference landmarks in image frame')
        plt.show()
    
    #--------------------------------------------------------------------------

    ### Cloud alignment via PCA

    # Prep
    pca = PCA()

    # Fit PCA model to data
    pca.fit(ref_lms)

    # Ensure that the sign of PCs is consistent with the image frame
    # Note: Given that the images are always acquired in the same orientations,
    #       a matching orientation can be ensured by finding the highest
    #       contributing image axis for each PC, and invert the PC if that
    #       contribution is negative. In other words, one ensures for each PC
    #       that the highest-contributing image axis is positively correlated 
    #       with the PC.
    
    # Find highest contributions of image axes to each PC
    # Note: This asks "which image axis contributes the most to this PC?"
    max_weights = np.argmax(np.abs(pca.components_),axis=1)
    
    # Get the signs of the highest contributions
    signs = np.sign(pca.components_[np.arange(pca.components_.shape[0]),max_weights])
    
    # Using the signs, flip those PCs where the sign is negative
    pca.components_ = pca.components_ * signs[:, np.newaxis]

    # Match the order of PCs to the order of image dimensions (zyx)
    # Note: Following the transform, the PCs will be sorted according to
    #       explained variance. Instead, they should be sorted in order of the
    #       highest contributing image dimension.
    
    # Find indices for zyx-sorting of transformed data   
    # Note: This asks "which PC is most contributed to by this image axis?"
    zyx_sort = np.argmax(np.abs(pca.components_),axis=0)

    # Transform landmarks, sort according to zyx
    lms_tf     = pca.transform(lms)[:,zyx_sort]
    ref_lms_tf = pca.transform(ref_lms)[:,zyx_sort]

    # Get PCs and explained variance to report
    PCs = np.copy(pca.components_.T)
    PCvars = np.copy(pca.explained_variance_ratio_)

    # Print results
    if verbose:
        print '\n  PCs:'
        print   '   ', str(PCs).replace('\n','\n    ')
        print   '  Explained variance:'
        print   '   ', str(PCvars)

    # Plot results
    if show:
        plt.scatter(lms_tf[:,2], lms_tf[:,1], c=lms_tf[:,0], 
                    cmap='viridis')
        plt.title('Channel landmarks in matched frame')
        plt.show()
        plt.scatter(ref_lms_tf[:,2], ref_lms_tf[:,1], c=ref_lms_tf[:,0], 
                    cmap='viridis')
        plt.title('Reference landmarks in matched frame')
        plt.show()
        
    
    #--------------------------------------------------------------------------

    ### Normalize z
    ### ...by scaling the 1st and 99th percentile to 0 and 1, respectively.

    # Get percentiles
    mem_bot = np.percentile(ref_lms_tf[:,0],1)
    mem_top = np.percentile(ref_lms_tf[:,0],99)

    # Scale
    lms_tf[:,0]     = (lms_tf[:,0]     - mem_bot) / (mem_top - mem_bot)
    ref_lms_tf[:,0] = (ref_lms_tf[:,0] - mem_bot) / (mem_top - mem_bot)        


    #--------------------------------------------------------------------------
    
    ### Additional Measure: Distance from Lumen
    
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
    
    # Change from pixels to resolution
    lumen = lumen * np.array(res)
    
    # Get Euclidean distance from lumen
    lum_dist_lms = np.sqrt(np.sum((lms-lumen)**2.0, axis=1))

    # Transform to PCA space
    lumen_tf = pca.transform(lumen.reshape(1,-1))[:,zyx_sort].squeeze()

    # Normalization of z
    lumen_tf[0] = (lumen_tf[0] - mem_bot) / (mem_top - mem_bot)

    # Report
    if verbose:
        print '  Lumen (raw & tf):'
        print '   ', lumen
        print '   ', lumen_tf

    # Plot to double-check
    if show:
        plt.scatter(ref_lms[:,2], ref_lms[:,1], c=ref_lms[:,0], 
                    cmap='viridis')
        plt.scatter(lumen[2], lumen[1], c='r', s=100)
        plt.title('Reference landmarks in image frame (with lumen)')
        plt.show()
        plt.scatter(ref_lms_tf[:,2], ref_lms_tf[:,1], c=ref_lms_tf[:,0], 
                    cmap='viridis')
        plt.scatter(lumen_tf[2], lumen_tf[1], c='r', s=100)
        plt.title('Reference landmarks in matched frame (with lumen)')
        plt.show()
        

    #--------------------------------------------------------------------------

    ### Return results
    
    return lms, lms_tf, lum_dist_lms


#------------------------------------------------------------------------------



