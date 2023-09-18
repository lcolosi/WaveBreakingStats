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

