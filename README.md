## Analysis and data repository for Terada et al., Nature (2021)
This repository contains code and links to dataset from Terada et al., Nature (2021). This code has been refactored to improve reproducibility.

### Prerequisities
* [MATLAB](https://ch.mathworks.com/products/matlab.html) (R2019a)
* [Signal Processing Toolbox](https://ch.mathworks.com/products/signal.html)
* [Statistics and Machine Learning Toolbox](https://ch.mathworks.com/products/statistics.html)

### Description

#### Detection_HSEs.m:
HSE_raw = Detection_HSEs(__) returns cell array to contain data matrices. Row of matrix correspond to individual region of intrests (ROIs) and columns correspond to time frames from event onset. Each column contains dF/F of individual ROIs at each timebin. The length of cell array corresponds the number of detected events in dataset.

* 'PL_raw' - matrix (ROI-by-timebin) of dF/F.
* 'mvFrame' - length of frames for moving avarage to normalize dF/Fs.
* 'position' - matrix containg position and velocity at each time frame.
* 'tWin' - number of frames to additionally remove before/after running events. (Any frames during running events are removed.)
* 'Ripple_onsets' - time frames of each ripple initiation
* 'latency_threshold' - criterion for ripple-associated or non-associated HSEs
* 'binWin' - number of frames to additionally extract around HSEs


#### TriggerDetect.m:
Crossings = TriggerDetect(__) returns cell array to contain two vectors and a matrix of them. A vector contains time frames of event initiation and another contains time frame of event terminations.

* 'Signal' - vector of data containing prominent events like spikes 
* 'UpThresh' - threshold for upward direction
* 'DownThresh' - threshold for downward direction


#### pairwise_reactivation.m
ccgZlagN = pirwise_reactivation(__) returns cross-correlation at zero-lag as the vector. The length corresponds position bins. 

* 'HSE_raw' - cell array detected by 'Detection_HSEs.m'
* 'HSE_coeff' - correlation coefficients between HSE (average df/f) and each ROI
* 'th_coeff' - threshold for correlation coefficient
* 'myPSTH' - place fields of each ROI
* 'beltLength' - length of a treadmill belt
* 'binsize' - binsize for smoothing
* 'binshift' - binsize to shift smoothing window

#### Authors
*st3166 AT columbia DOT edu

#### Acess to the dateset:

The dataset will be made available soon
