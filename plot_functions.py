#!/usr/bin/env python

""" plot_functions.py: suite of visualizations and kinematic studies for TrackML events """ 

__author__  = "Gage DeZoort"
__version__ = "1.0.0"
__status__  = "Development"

import os
import math
import trackml
import numpy as np
import pandas as pd
from cycler import cycler
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import mpl_toolkits.mplot3d.art3d as art3d
from trackml.dataset import load_event
from trackml.dataset import load_dataset


def plotSingleHist(data, x_label, y_label, bins, weights=None, title='', color='blue'):
    """ plotSingleHist(): generic function for histogramming a data array
    """
    bin_heights, bin_borders, _ = plt.hist(data, bins=bins, color=color, weights=weights)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders)/2.
    print bin_heights
    print bin_centers
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.title(title,    fontsize=16)
    plt.show()

def plotBinnedHist(x, y, x_label, y_label, nbins, title='', color='blue'):
    """ plotBinnedHist(): generic function for plotting a binned y vs. x scatter plot
    					  similar to a TProfile
    """
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)/np.sqrt(n)
    plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, color=color, marker='.', ls='')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()

def plotHeatMap(x, y, x_label, y_label, x_bins, y_bins, weights=None, title=''):
    """ plotHeatMap(): generic function for plotting a heatmap of y vs. x with 
                       pre-specified weights
    """
    plt.hist2d(x, y, bins=(x_bins,y_bins), weights=weights, cmap=plt.cm.jet)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Normalized $P_T$')
    plt.title(title)
    plt.show()

def plotTrack(track):
	""" plotTrack(): plot a track in a simple 3D grid 
	"""
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot(track['tx'], track['ty'], track['tz'], lw=0.5, c='skyblue')
        ax.scatter3D(track['tx'], track['ty'], track['tz'], 
                     c=track['tR'], cmap ='viridis', marker='h', s=30)
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')
        ax.set_zlabel('z [mm]')
        plt.show()

def plotTrackOverLayers(track, volume_ids):
    """ plotTrackOverLayers(): plot a track and the detector layers
                               it hits
    """

    detectors = pd.read_csv('../detectors.csv')
    detectors['xyz'] = detectors[['cx', 'cy', 'cz']].values.tolist()
    
    volumes = detectors.groupby('volume_id')['xyz'].apply(list).to_frame()	
    accept_volumes = detectors[detectors.volume_id.isin(volume_ids)]
    
    x_min, x_max = accept_volumes['cx'].min(), accept_volumes['cx'].max()
    y_min, y_max = accept_volumes['cy'].min(), accept_volumes['cy'].max()
    z_min, z_max = accept_volumes['cz'].min(), accept_volumes['cz'].max()
    
    volumes_layers = accept_volumes.groupby(['volume_id','layer_id'])['xyz'].apply(list).to_frame()
    fig = plt.figure(figsize=plt.figaspect(0.9))
    
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    #ax.set_prop_cycle(cycler('color', ['forestgreen', 'firebrick', 'royalblue', 'indigo']))
    
    ax.plot(track['tx'], track['ty'], track['tz'], lw=0.5, c='skyblue')
    ax.scatter3D(track['tx'], track['ty'], track['tz'],
                 c=track['tR'], cmap='viridis', marker='h', s=30)
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    
    num_regions = volumes_layers.shape[0]
    for (i, row) in volumes_layers.iloc[:num_regions+1].iterrows():
        xyz = np.array(row['xyz'])
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        ax.plot(x,y,z, linestyle='', marker='h', markersize='0.5', color='mediumslateblue', alpha=0.5)
        ax.text(x[0], y[0], z[0], str(i), None, size=5)
    plt.show()
        







