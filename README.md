## GNN Algorithms for Tracklet Finding in the CMS/ATLAS Pixel Detectors 
### Introduction 
This code uses data from the CERN [TrackML Challenge](https://www.kaggle.com/c/trackml-particle-identification/overview) which can be downloaded from Kaggle. This repository includes the python library to easily read the data into a dataframe in the `trackml-library` folder. Upon first downloading this repository you must install this library by running 
``` 
pip install /path/to/repository
```
### Data
TrackML events are separated into four files: cells, hits, particles, and truth. The **data** folder contains event0000010000 from the TrackML dataset. 

### Visualization and Kinematic Studies
Several visualization scripts and kinematic studies are available in the **visualization_scripts** folder. 
* **skim_data.py**: produce skim event files (keeps only hits/truth in the pixel detector, removes all cell information, and adds detector/cylindrical coordinates to remaining truth hits) 
* **plot_functions.py**: contains generic functions for plotting histograms (block, binned), scatterplots (errorbars, overlayed plots), heat maps, tracks in 3D space, tracks overlapped with the modules they hit, and entire regions of the detector.
* **segment_detector.py**: returns the detector dataframe where each hit has been sorted into one of N phi bins
* **analyze_tracks.py**: analyze track quality vs. non-quality, pt per track distributions, number of hits per track distributions, number of layers hit per track distributions, track length distributions, dEta/dPhi/dR per track distributions

### Data Structures
The goal of pre-processing the data is to create graphs, which are namedtuples of matrices X, Ri, Ro, and y. X is the feature vector, which contains the cylindrical position (r, phi, z) of each hit. Ri and Ro are segment matrices, each of which have nHits rows and nSegments columns. Element Ri_{hs} of Ri is 1 if segment s is incoming to hit h, and 0 otherwise. Likewise, element Ro_{hs} of Ro is 1 if segment s is outgoing from hit h, and 0 otherwise. y is the segment truth vector, which is a vector of length nSegments containing 0 entries for false segments and 1 entries for true segments. The **graph.py** file defines graphs and some corresponding loading/saving functions. 

### Measurements
We have defined several graph construction performance metrics, which are implemented in this folder. Among them are the *segment efficiency*, which is defined as sum(y)/len(y), and *truth efficiency*, which is the number of true segments selected divided by the total number of true segments contained in the dataset. In the latter case, it is necessary to calculate the correct number of truth segments for each pre-processing strategy. For example, in the layer pairs pre-processing scheme (in which only one hit per layer per particle is kept), the true number of hits per particle is simply nLayersHit-1. 

* **graph_efficiency.py**: calculate and compares segment efficiency and truth efficiency across two different pre-processing strategies for different pt cuts 
* **truth/generate_truth_LP.py**: calculates the true number of segments that should be caputured by the layer pair strategy (for a range of pt cuts)
* **truth/generate_truth_LPP.py**: calculates the true number of segments that should be captured by the layer pair+ strategy (for a range of pt cuts)

