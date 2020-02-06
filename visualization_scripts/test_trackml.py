import os
import math
import trackml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trackml.dataset import load_event
from trackml.dataset import load_dataset

test_evt = 'event000001000'
hits, cells, particles, truth = load_event(os.path.join('train_100_events', test_evt))

# Add pt, eta, phi coordinates 
particles['p']   = np.sqrt(particles['px']**2 + particles['py']**2 
                           + particles['pz']**2)
particles['pt']  = np.sqrt(particles['px']**2 + particles['py']**2)
particles['eta'] = np.arctanh(particles['pz']/particles['p'])
particles['phi'] = np.arctan2(particles['py'],particles['px']) # range [-pi,pi] 
particles.loc[particles.phi < 0, 'phi'] += 2*np.pi # range [0,2pi] w.r.t. x-axis

# Grab different columns to plot
pt_large  = particles[particles.pt >= 5].pt
pt_small  = particles[particles.pt < 5].pt
eta = np.array(particles.eta)
phi = np.array(particles.phi)
pt_eta_phi = particles[['pt','eta','phi']]
pt_sm_eta_phi = pt_eta_phi[pt_eta_phi.pt < 5]
pt_lg_eta_phi = pt_eta_phi[pt_eta_phi.pt >= 5]
pt_large_bins = np.linspace(5, 105, num=50)
pt_small_bins = np.linspace(0, 5,   num=50) 
phi_bins = np.linspace(0, 2*np.pi, num=15)
eta_bins = np.linspace(-5, 5, num=15)

def plotSingleHist(data, x_label, y_label, bins, weights=None, density=False,
                   title='', color='blue'):
    bin_heights, bin_borders, _ = plt.hist(data, bins=bins, color=color, weights=weights)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders)/2.
    plt.legend(loc = 'upper right')
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.title(title,    fontsize=16)
    plt.show()

def plotHeatMap(x, y, x_label, y_label, x_bins, y_bins, weights=None, title=''):
    # x_heights, x_edges = np.histogram(x, x_bins, weights=weights)
    # y_heights, y_edges = np.histogram(y, y_bins, weights=weights)
    plt.hist2d(x, y, bins=(x_bins,y_bins), weights=weights, cmap=plt.cm.jet)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Normalized $P_T$')
    plt.title(title)
    plt.show()


print ("Question: where does all the pT go?" )
print ("... plotting pt_small histogram")
plotSingleHist(pt_small, "$p_T$ [GeV]", "Counts", pt_small_bins, title="$p_T<5$ GeV", color='mediumslateblue')

print ("... plotting pt_small histogram")
plotSingleHist(pt_large, "$p_T$ [GeV]", "Counts", pt_large_bins, title="$p_T\geq 5$ GeV", color='firebrick')

norm = float(pt_sm_eta_phi['pt'].sum())
weights = np.array(pt_sm_eta_phi['pt'].div(norm))

print( "... plotting pt_small vs. phi histogram")
plotSingleHist(pt_sm_eta_phi['phi'], '', '', bins=25, weights=weights)

print( "... plotting pt_small vs. eta histogram")
plotSingleHist(pt_sm_eta_phi['eta'], '', '', bins=25, weights=weights)

print ("... plotting pt vs. phi,eta heatmap")
plotHeatMap(np.array(pt_sm_eta_phi['phi']), np.array(pt_sm_eta_phi['eta']), "$\phi$", "$\eta$",
            x_bins=phi_bins, y_bins=eta_bins, weights=weights)


