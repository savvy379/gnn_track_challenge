#!/usr/bin/env python

import os
import sys
sys.path.append("../")

import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib.pyplot as plt
from matplotlib import rc

from data_structures.graph import Graph, load_graph
import visualization_scripts.plot_functions as pf

# hard-coded parameters
n_events = 100
pt_cuts  = ['0p5', '0p6', '0p75', '1', '1p5', '2', '2p5', '3', '4', '5']

# configure matplotlib
font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 12}
rc('font', **font)
rc('text', usetex=True)

def get_graphs(d):
    """ get all graphs from a directory
          return [("event<id1>", graph1, size), ...]
    """
    files = os.listdir(d)
    return [(f.split('_')[0], load_graph(d+f), os.path.getsize(d+f)/10**6) 
            for f in files]

def get_segment_efficiency(g):
    """ segment efficiency = sum(y)/len(y)
    """
    return [np.sum(g[i][1].y)/(g[i][1].y).shape[0]
            for i in range(n_events)]

def get_shape(g):
    """ graph shape = (n_nodes, n_edges)
    """
    print(np.array(g).shape)
    n_nodes = [g[i][1].X.shape[0] for i in range(n_events)]
    n_edges = [g[i][1].Ri.shape[1] for i in range(n_events)]
    return np.array([n_nodes, n_edges])

def get_truth_efficiency(i, gr, tr):
    """ truth efficiency = sum(y)/n_total_true
    """
    truth_effs = [np.sum(g[1].y)/float(tr.loc[tr['evt_id']==g[0], pt_cuts[i]]) 
                  for g in gr]
    return truth_effs

def print_summary(eff, nodes, edges, tag=""):
    print(tag)
    print("  ==> Average n_nodes = ", np.round(nodes[0], decimals=2),
          "+/-", np.round(nodes[1], decimals=2))
    print("  ==> Average n_edges = ", np.round(edges[0], decimals=2),
          "+/-", np.round(edges[1], decimals=2))
    print("  ==> Segment Efficiency", np.round(eff[0], decimals=3),
          "+/-", np.round(eff[1], decimals=3))    

def read_truth(file_name):
    """ truth file contains the true number of segments per
        pt cut for each file
    """
    truth_file = open(file_name, 'r')
    lines = [i.split(" ") for i in truth_file]
    truth_info = df({'evt_id' : [], '0p5' : [], '0p6' :[],
                     '0p75' : [], '1' : [], '1p5' : [], 
                     '2' : [], '2p5' : [], '3' : [], 
                     '4' : [], '5' : []})
    for line in lines:
        truth_info = truth_info.append({'evt_id' : line[0],
                                        '0p5'    : int(line[1]),
                                        '0p6'    : int(line[2]),
                                        '0p75'   : int(line[3]),
                                        '1'      : int(line[4]),
                                        '1p5'    : int(line[5]),
                                        '2'      : int(line[6]),
                                        '2p5'    : int(line[7]),
                                        '3'      : int(line[8]),
                                        '4'      : int(line[9]),
                                        '5'      : int(line[10])},
                                       ignore_index=True)
        
    return truth_info


truth_info = [read_truth("truth/truth_LP.txt"), read_truth("truth/truth_LPP.txt")]
data0_dirs  = ['/tigress/jdezoort/prep_100events_' + pt + '/' 
              for pt in pt_cuts]
data1_dirs = ['/tigress/sthais/prep_graphs/layer_pair_plus/lpp_' + pt + '/'
              for pt in pt_cuts]

dirs = np.array([data0_dirs, data1_dirs])

to_hist_0 = df({'seg_eff'   : [], 'seg_eff_er'   : [],
                'truth_eff' : [], 'truth_eff_er' : [],
                'n_segs'    : [], 'n_segs_er'    : [],
                'n_nodes'   : [], 'n_nodes_er'   : [],
                'size'      : [],  'size_er'     : []})

to_hist_1 = df({'seg_eff'   : [], 'seg_eff_er'   : [],
                'truth_eff' : [], 'truth_eff_er' : [],
                'n_segs'    : [], 'n_segs_er'    : [],
                'n_nodes'   : [], 'n_nodes_er'   : [],
                'size'      : [],  'size_er'     : []})

for d in range(dirs.shape[0]):

    print("Examining directory", dirs[d])

    for i in range(dirs.shape[1]):

        data_graphs = get_graphs(dirs[d][i])
        data_shapes = get_shape(data_graphs)
        
        
        truth_effs = get_truth_efficiency(i, data_graphs, truth_info[d])
        seg_effs   = get_segment_efficiency(data_graphs)
        
        avg_seg_eff   = [np.mean(seg_effs),       np.sqrt(np.var(seg_effs))]
        avg_truth_eff = [np.mean(truth_effs),     np.sqrt(np.var(truth_effs))]
    
        avg_nodes     = [np.mean(data_shapes[0]), np.sqrt(np.var(data_shapes[0]))]
        avg_edges     = [np.mean(data_shapes[1]), np.sqrt(np.var(data_shapes[1]))]
        avg_size      = [np.mean([g[2] for g in data_graphs]),
                         np.sqrt(np.var([g[2] for g in data_graphs]))]

        data_tag = " ***** pt=" + pt_cuts[i] + " data ***** "
        print_summary(avg_seg_eff, avg_nodes, avg_edges, tag=data_tag)

        if (d==0):
            to_hist_0 = to_hist_0.append({'seg_eff'      : avg_seg_eff[0],
                                          'seg_eff_er'   : avg_seg_eff[1],
                                          'truth_eff'    : avg_truth_eff[0],
                                          'truth_eff_er' : avg_truth_eff[1],
                                          'n_segs'       : avg_edges[0],
                                          'n_segs_er'    : avg_edges[1],
                                          'n_nodes'      : avg_nodes[0],
                                          'n_nodes_er'   : avg_nodes[1],
                                          'size'         : avg_size[0],
                                          'size_er'      : avg_size[1]},
                                         ignore_index=True)
        if (d==1):
            to_hist_1 = to_hist_1.append({'seg_eff'      : avg_seg_eff[0],
                                          'seg_eff_er'   : avg_seg_eff[1],
                                          'truth_eff'    : avg_truth_eff[0],
                                          'truth_eff_er' : avg_truth_eff[1],
                                          'n_segs'       : avg_edges[0],
                                          'n_segs_er'    : avg_edges[1],
                                          'n_nodes'      : avg_nodes[0],
                                          'n_nodes_er'   : avg_nodes[1],
                                          'size'         : avg_size[0],
                                        'size_er'      : avg_size[1]},
                                         ignore_index=True)
            
pt_cuts = np.array([0.5, 0.6, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5])
pf.plotXYXY(pt_cuts, np.array(to_hist_0['seg_eff']), 'LP', 
         pt_cuts, np.array(to_hist_1['seg_eff']), 'LP+',
         '$p_{T}$ Cut [GeV]', 'True/Total Selected Segments', 
         yerr1=np.array(to_hist_0['seg_eff_er']),
         yerr2=np.array(to_hist_1['seg_eff_er']),
         title='segment_efficiency_comparison',
         color1='mediumslateblue', color2='mediumseagreen')

pf.plotXYXY(pt_cuts, np.array(to_hist_0['truth_eff']), 'LP',
          pt_cuts, np.array(to_hist_1['truth_eff']), 'LP+',
         '$p_{T}$ Cut [GeV]', '(Selected True)/(Total True) Segments',
          yerr1=np.array(to_hist_0['truth_eff_er']),
          yerr2=np.array(to_hist_1['truth_eff_er']),
          title='truth_efficiency_comparison',
          color1='mediumslateblue', color2='mediumseagreen')

pf.plotXY(pt_cuts, np.array(to_hist_0['seg_eff']), '$p_{T}$ Cut [GeV]',
       'True/Total Selected Segments', yerr=np.array(to_hist_0['seg_eff_er']), 
       color='mediumslateblue', title='segment_efficiency_LP')

pf.plotXY(pt_cuts, np.array(to_hist_0['truth_eff']), '$p_{T}$ Cut [GeV]',
       '(Selected True)/(Total True) Segments', yerr=np.array(to_hist_0['truth_eff_er']), 
       color='mediumslateblue', title='truth_efficiency_LP')

pf.plotXY(pt_cuts, np.array(to_hist_1['seg_eff']), '$p_{T}$ Cut [GeV]',
       'True/Total Selected Segments', yerr=np.array(to_hist_1['seg_eff_er']), 
       color='mediumseagreen', title='segment_efficiency_LPP')

pf.plotXY(pt_cuts, np.array(to_hist_1['truth_eff']), '$p_{T}$ Cut [GeV]',
       '(Selected True)/(Total True) Segments', yerr=np.array(to_hist_1['truth_eff_er']), 
       color='mediumseagreen', title='truth_efficiency_LPP')



"""
plotXY(pt_cuts, np.array(build_time)/float(n_events),
       '$p_{T}$ Cut [GeV]', 'Avg. Construction Time [s]', color='indigo',
       title='construction_time')
plotXY(pt_cuts[4:], np.array(build_time)[4:]/float(n_events),
       '$p_{T}$ Cut [GeV]', 'Avg. Construction Time [s]', color='indigo',
       title='construction_time_zoomed')
plotXY(pt_cuts, np.array(to_hist['size']), '$p_{T}$ Cut [GeV]', 
       'Size [MB]', yerr=np.array(to_hist['size_er']), color='indigo',
       title='size')
plotXY(pt_cuts[4:], np.array(to_hist['size'])[4:], '$p_{T}$ Cut [GeV]',
       'Size [MB]', yerr=np.array(to_hist['size_er'])[4:], color='indigo',
       title='size_zoomed')
plotXY(pt_cuts, np.array(to_hist['seg_eff']), '$p_{T}$ Cut [GeV]', 
       'True/Total Selected Segments', yerr=np.array(to_hist['seg_eff_er']), color='dodgerblue',
       title='segment_efficiency')
plotXY(pt_cuts, np.array(to_hist['truth_eff']), '$p_{T}$ Cut [GeV]', 
       'Selected/Total True Segments', yerr=np.array(to_hist['truth_eff_er']), color='dodgerblue',
       title='truth_efficiency')
plotXY(pt_cuts, np.array(to_hist['n_segs']), '$p_{T}$ Cut [GeV]', 
       'Segments per Graph', yerr=np.array(to_hist['n_segs_er']), color='seagreen',
       title='n_segments')
plotXY(pt_cuts, np.array(to_hist['n_nodes']), '$p_{T}$ Cut [GeV]',
       'Nodes per Graph', yerr=np.array(to_hist['n_nodes_er']), color='seagreen',
       title='n_nodes')
"""
