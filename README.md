# Wave Breaking Statistics Code Repository
Luke Colosi | lcolosi@ucsd.edu 

# Overview

Here, we provide the MATLAB code for computing wave breaking statistics based on the seminal work of Phillips 1986.
The MATLAB scripts take georeferenced video camera images of the ocean surface and estimate the lambda of c omnidirectional and directional distributions along with the moments of the distribution.

Video camera is apart of the Modular Aerial Sensing System (MASS), a compact airborne remote sensing instrument package developed by the Air-Sea Interactions Lab at SIO for observing surface and upper ocean processes (see Melville et al. 2016 for more information).
More specifically, the video camera is an IO Industries FLARE 12M125-CL camera with 4096 by 3072 pixel resolution, 10 bit, and  5-Hz sampling rate.     

# Outline of file system structure 

Each experiment that the MASS is deployed has the following file system structure (here I use S-MODE IOP1 as the example):

```
SMODE_2022/
├── MAT
├── PROCESSED
│   ├── DoppVis
│   ├── HYPER_Proc
│   ├── KT19
│   ├── LIDAR
│   ├── LIDAR_SDC
│   ├── Trajectory
│   └── VIDEO
│       ├── 16bit_TIF_Frames
│       ├── IE_Camera_Event_EO
│       └── Trimble
├── PROGRAMS
└── RAW
```

`SMODE_2022` is the root of this experiment's file tree and contains the `MAT`, `PROCESSED`, `PROGRAMS`, and `RAW` subdirectories along with others depending on the experiment.
The most important subdirectories are the following: 

`Processed`: Contains the processed data from each of the MASS instruments. 
`PROGRAMS`: Contains the MATLAB code (main scripts and functions) for processing and analyzing MASS data. 
`RAW`: Contains the raw (unprocessed) data from each of the MASS instruments. 

Focusing on just the video camera, under `SMODE_2022/PROCESSED/VIDEO/`, we have the following subdirectories: 

`16bit_TIF_FRAMES`: Contains raw (nongeoreferenced) video images (found in images/ subdirectory). 
`IE_Camera_Event_EO`: 
`Trimble`: Contains georeferenced video images. 

Each of these subdirectory contains the directories or files that are named based on the data of each flight.
Right now, all level-2 data (derived quantities from georeferenced vide images) as well as quality control figures and a bunch of development files are saved inside `SMODE_2022/PROCESSED/VIDEO/Trimble/date_of_flight/` directories. 
This could be cleaned up by creating another subdirectory in `SMODE_2022/PROCESSED/VIDEO/` dedicated to level-2 data and the quality control figures illustrating the processing steps.    

# 


