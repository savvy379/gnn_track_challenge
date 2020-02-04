#!/usr/bin/env python

""" skim_data.py: produce skimmed datasets""" 

__author__  = "Gage DeZoort"
__version__ = "1.0.0"
__status__  = "Development"

import os
import math
import trackml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d

from trackml.dataset import load_event
from trackml.dataset import load_dataset
from os import listdir
from os.path import isfile, join

def getArgs():
	""" getArgs(): options -N <NumberOfFiles> and -d <DirectoryWithEvents>
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-input", default="TrainEvents", type=str, help="Directory containing event files.")
	parser.add_argument("-output", default="TrainEvents", type=str, help="Directory to write skimmed files.")
	return parser.parse_args()

def addDetectorCoords(df):
	""" addDetectorCoords(): input a dataframe and return it with 
                             additional pt, eta, and phi coords
	"""
	df['p']     = np.sqrt(df['px']**2 + df['py']**2 + df['pz']**2)
	df['pt']    = np.sqrt(df['px']**2 + df['py']**2)
	df['eta']   = np.arctanh(df['pz']/df['p'])
	df['phi']   = np.arctan2(df['py'], df['px']) # range [-pi, pi]
	df.loc[df.phi < 0, 'phi'] += 2*np.pi # range [0, 2pi]
	return df


args = getArgs()
path = args.input
out_dir = args.output
if out_dir[-1] == '/': 
	out_dir = out_dir[:-1]
if path[-1] == '/':
	path = path[:-1]
if not os.path.exists(out_dir):
	print "mkdir", out_dir
	os.makedirs(out_dir)
print "Writing to", out_dir

# transform filenames into event_strings
files = [f for f in listdir(path) if isfile(join(path, f))]
files = [f.split('.')[0] for f in files if "truth" in f]
files = [f.split('-')[0] for f in files]

num_processed = 0
for f in files:
	
	print " -> Processing", path + "/" + f + "*"
	hits, cells, particles, truth = load_event(os.path.join(path, f))
	
	# drop hits outside of the pixel detector 
	hits = hits[(hits.volume_id < 10)]
	hit_ids = hits.hit_id.unique()
	truth = truth[truth.hit_id.isin(hit_ids)]

	# add particle physics coordinates 
	truth['tR'] = np.sqrt(truth['tx']**2 + truth['ty']**2 + truth['tz']**2)
	truth = truth.rename(columns={'tpx':'px', 'tpy':'py', 'tpz':'pz'})
	truth = addDetectorCoords(truth)
	particles = addDetectorCoords(particles)

	# write skimmed dataframes to files
	truth.to_csv(out_dir+"/"+f+"-truth.csv")
	hits.to_csv(out_dir+"/"+f+"-hits.csv")
	particles.to_csv(out_dir+"/"+f+"-particles.csv")
	
	# note that this cell file is a dummy
	cells[cells.hit_id==1].to_csv(out_dir+"/"+f+"-cells.csv")
	num_processed += 1

print "Completed", num_processed, "event skims!"






