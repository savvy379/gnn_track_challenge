#!/usr/bin/env python

""" segment_detector.py: functions for binning the pixel detector """

__author__  = "Gage DeZoort"
__version__ = "1.0.0"
__status__  = "Development"

import math
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

def segment_detector(detectors, N):
  ''' segment_detector(): analyze detector segmenting, return detector with 
                          phi coordinates and phi bins '''

  unique_volumes = [7,8,9]
  pixel_detector = detectors[detectors.volume_id.isin(unique_volumes)]

  # add phi to pixel detector module coordinates
  pixel_detector['phi'] = np.arctan2(pixel_detector['cy'], pixel_detector['cx'])
  pixel_detector.loc[pixel_detector.phi < 0, 'phi'] += 2*np.pi
  
  # create phi bins for the modules
  n_bins = N
  phi_bins = np.linspace(0, 2*np.pi, n_bins+1)
  total_modules_per_bin = np.zeros(n_bins)
  bins = np.digitize(pixel_detector['phi'], phi_bins)
  pixel_detector['bin'] = bins
  print pixel_detector['bin'].unique()
  print pixel_detector.head(500)

  print "\n*** Beginning Segmentation Analysis *** \n"
  print "    --> Number of phi bins:", n_bins
  print "    --> Specific phi_bins/pi:\n       ", phi_bins/np.pi
  
  # loop over each pixel detector volume
  for volume in unique_volumes:
    modules = pixel_detector[pixel_detector.volume_id == volume]
    unique_layers = modules.layer_id.unique()
    print "    --> Unique layers in Vol", str(volume) + ":", unique_layers 
    
    # number modules per layer
    module_counts = modules.groupby('layer_id').count()['module_id']
    
    for layer in unique_layers:
      
      phi_coords = np.array(modules[modules.layer_id == layer]['phi'])
      phi_binned_modules = np.histogram(phi_coords, bins=phi_bins)[0]
      total_modules_per_bin = total_modules_per_bin + phi_binned_modules
      
  print "    --> Total modules per phi bin:\n       ", total_modules_per_bin
  pixel_detector.reset_index(drop=True)
  return pixel_detector
