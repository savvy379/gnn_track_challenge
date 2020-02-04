#!/usr/bin/env python

""" analyze_tracks.py: perform a basic analysis of track features """

__author__  = "Gage DeZoort"
__version__ = "1.0.0"
__status__  = "Development"

import os
import sys
from os import listdir
from os.path import isfile, join
sys.path.append("..")
import math
import trackml
import numpy as np
import pandas as pd
import plot_functions as pf
from trackml.dataset import load_event
from trackml.dataset import load_dataset

def forwardProgress(df):
  """ forwardProgress(): check that track makes radial progress
  """
  radial_steps = df['tR'].diff().dropna()
  positiveStep = False
  negativeStep = False
  for step in radial_steps:
    if abs(step) > 10:
      if step < 0: negativeStep = True
      if step > 0: positiveStep = True
      if negativeStep and positiveStep: return False
  return True

N_to_analyze = 100    # total number of events to analyze
N_total = 0              # number of analyzed events
N_good = 0               # number of good events
N_backtracks = 0         # number of backtracks
N_sparse = 0             # number of events with nHits < 3

to_hist = pd.DataFrame({'pt' : [], 'dPhi' : [], 'dEta' : [], 
			'dR' : [], 'len' : [], 'nHits' : []})

files = [f for f in listdir("skimmed_events/") if isfile(join("skimmed_events/", f))]
files = [f.split('.')[0] for f in files if "truth" in f]
files = [f.split('-')[0] for f in files]

for evt_num in files:

  print "Processing", evt_num
  hits, cells, particles, truth = load_event(os.path.join('skimmed_events', evt_num))

  particle_ids = truth.particle_id.unique()

  for id in particle_ids:
		
    if int(id) == 0: continue
    particle = particles[particles.particle_id == id]
    
    # basic pT cut
    if (particle.iloc[0,:]['pt'] < 5): continue
    
    truth_track = truth[truth.particle_id == id]
    layers_hit  = hits.loc[hits['hit_id'].isin(truth_track['hit_id'])]
		
    #pf.plotTrackOverLayers(truth_track, np.array([7,8,9]))
    N_total += 1

    skipEvent = False
    nHits = truth_track.shape[0]
    if nHits < 3:
      N_sparse += 1
      skipEvent = True

    if not forwardProgress(truth_track):
      N_backtracks += 1
      skipEvent = True

    if skipEvent: continue

    first_hit = truth_track.iloc[0,:]
    last_hit = truth_track.iloc[-1,:]
    pt = particle.iloc[0,:]['pt']
    dEta = abs(last_hit['eta'] - first_hit['eta'])
    dPhi = min(2*np.pi-abs(last_hit['phi'] - first_hit['phi']),
	       abs(last_hit['phi'] - first_hit['phi']))
    dR = np.sqrt(dEta**2 + dPhi**2)
    len  = np.sqrt((last_hit['tx'] - first_hit['tx'])**2 + 
		   (last_hit['ty'] - first_hit['ty'])**2 + 
		   (last_hit['tz'] - first_hit['tz'])**2)
    
    
    to_hist = to_hist.append({'pt' : pt, 'dEta' : dEta, 'dPhi': dPhi,
			      'dR' : dR, 'len' : len,
			      'nHits' : nHits}, ignore_index=True)
    
    N_good += 1
    if N_total >= N_to_analyze: break

  if N_total >= N_to_analyze: break

print " --> There were", N_backtracks, "/", N_total, "backtracks!"
print " --> There were", N_sparse, "/", N_total, "events with <3 hits!"
print " --> Good events analyzed:", N_good, "/", N_total

#pf.plotSingleHist(np.array(to_hist['pt']), '$p_T$ [GeV]', 'Counts', bins=[0,1,2,5,10,20,40,60,100,150], color='mediumseagreen')
pf.plotSingleHist(np.array(to_hist['pt']), '$p_T$ [GeV]', 'Counts', bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5], color='mediumseagreen')
pf.plotSingleHist(np.array(to_hist['nHits']), 'nHits', 'Counts', bins=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], color='sandybrown')
pf.plotSingleHist(np.array(to_hist['len']), 'Track Length [mm]', 'Counts', bins=15, color='darkviolet')
pf.plotBinnedHist(np.array(to_hist['pt']), np.array(to_hist['nHits']), '$p_T$ [GeV]', 'nHits', 70, title='nHits vs. $p_T$ for all tracks', color='sandybrown')
pf.plotBinnedHist(np.array(to_hist['pt']), np.array(to_hist['len']), '$p_T$ [Gev]', 'Track Length [mm]', 70, color='darkviolet', title='Track Length vs. $p_T$ for all tracks')

pf.plotBinnedHist(np.array(to_hist['pt']), np.array(to_hist['dR']), '$p_T$ [GeV]', '$dR=\sqrt{d\eta^2+d\phi^2}$', 70, title='$dR$ vs. $p_T$ for all tracks', color='cornflowerblue')
pf.plotBinnedHist(np.array(to_hist['pt']), np.array(to_hist['dEta']), '$p_T$ [GeV]', '$d\eta$', 70, title='$d\eta$ vs. $p_T$ for all tracks', color='cornflowerblue')
pf.plotBinnedHist(np.array(to_hist['pt']), np.array(to_hist['dPhi']), '$p_T$ [GeV]', '$d\phi$', 70, title='$d\phi$ vs. $p_T$ for all tracks', color='cornflowerblue')




