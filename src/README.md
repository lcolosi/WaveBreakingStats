# Main.m Script
*Authors: Teodor Vrecica, Luc Lenain, and Luke Colosi*

Main script for estimating the $\Lambda(c)$ distribution and other wave-breaking statistics from visual imagery recorded by a video camera mounted onboard research aircraft. Images are georeferenced using a coupled GPS-IMU, which provides aircraft trajectory information.

This script imports the following data: 

1. **Raw (non-georeferenced) video camera images**: Consists of .tif files (organized by flight and track numbers) recorded by the IO Industries FLARE 12M125-CL camera. Images are collected at a 5-Hz sampling rate with a 4096 by 3072-pixel resolution. 

2. **Aircraft trajectory and attitude**: Consists of the EO text files that provide the position and attitude of the aircraft over time from the coupled GPS-IMU onboard the research aircraft. The coupled GPS-IMU consists of the Novatel SPAN LN200 Inertial Motion Unit (IMU) and a ProPak6 GPS receiver. This time series is interpolated to match the time steps of the video camera images.  

For specifics on video and trajectory data processing, see Nick Statom's documentation [here](https://docs.google.com/document/d/1qbaBH98IW1tJrMfxC6TKL-jQIcMPm7KK_KRrxk3rBb0/edit).

## Program structure 
Estimating $\Lambda(c)$ distribution includes the following steps: 

1. Identifying trajectory and visual imagery data for flight tracks.
2. Quality controlling data based on aircraft stability (i.e. how much the plane rolls, pitches, and changes heading).
3. Remove vignetting and computing parameters for identifying regions of high sun-glint.
4. Build and execute batch files for georeferencing processed images with the Trimble software. 
5. Determining the brightness threshold of breakers based on Kleiss and Melville (2011).  
6. Computing the $\Lambda(c)$ distribution and its moments for wave breaking statistics.

## Input parameters
### Select processes to run (0 or 1)
1. `option_plot` : Plot verification figures
2. `option_image_proc` : Process raw images and generate intermediate data products (unnecessary if these intermediate products are already generated and you just want to work lambda of c distribution code)
3. `option_globalOrlocal` : Check if a global brightness threshold is defined for the area (1), or a local one for each image (0).

### Date
1. `StartDate` : Start date of flight in UTC time (format: 'yyyymmdd')

### Directories
1. `dirRaw` : Path to the directory containing Raw (non-georeferenced) video images for a given flight. These are usually located in the `Images\` directory of the Trimble project. 
2. `dirProc` : Path to the directory containing the Trimble project directories and files. The Trimble project directories include `EO\`, `Images\`, `Output\`, `Project\`, and `Template\`. 
3. `dirV` : Path to the directory where quality control/verification plots are saved.
4. `dirOut` : Path to the directory where intermediate data products (e.g., .mat files) are saved.

### Flight Stability criteria parameters
1. `maxPer` : Maximum percent of flight track to be removed at the beginning and end of the track. For example, maxPer = [25,25] means that at a minimum, 50 percent (from 25% - 75%) of the flight track will be used in analysis with 0-25% and 75%-100% removed from the flight track.  
2. `sigRoll` : Maximum allowed roll standard deviations for a stable segment.    
3. `sigPitch` : Maximum allowed pitch standard deviations for a stable segment.
4. `sigHeading` : Maximum allowed heading standard deviations for a stable segment.
5. `Shift` : Shifts in flight track for finding a stable flight period
6. `Nstd` : Number of standard deviations of either roll, pitch, or heading that constitute an abrupt change in attitude of the plane
7. `tcheck` : Time interval between the jth roll/pitch/heading observation to check if there is an abrupt change in attitude shortly after the jth observation.

### Vignette removal parameters
1. `winSize` : Side length of the square kernel used in the 2D moving average for smoothing (low pass filtering) mean brightness of each pixel over the track period. 
2. `sigma_ff_m` : Standard deviation of the Gaussian smoothing filter for the 2D image flat-field correction used in calculating the mean and standard deviation pixel image brightness.
3. `sigma_ff_v` : Same as above, but for vignetting removal calculation.
4. `sigBrightness` : Maximum allowed standard deviation of meanOriginal to constitute constant brightness along the track.
5. `B_threshold` : Fraction of the pixels that are considered in the mean pixel brightness calculation.
6. `n_sigma` : Number of standard deviations above or below the median image brightness. The parameter is used for determining which pixels are considered for the calculation of the mean image brightness.
7. `stdMag` : Number of standard deviations away from the median value of the mean pixel brightness image. Used to determine the sun glint brightness threshold. 
   
### Trimble georeferencing project parameters
1. `camName` : Video Camera name 
2. `utmZone` : UTM zone for experiment
   
### Brightness threshold parameters
1. `localStep` :
2. `peakPercentage` :

## Tips for running the code 

1. This code can be run locally on your computer or remotely through whichever airseaserver the data is located on. Either way, your local computer needs to be connected to the UCSD VPN; use Cisco AnyConnect to do this. If running the code locally, connect to the remote server  through your file system explorer to gain remote access to the data (i.e., mount to the remote server on your local computer). If running the code on the airseaserver, use a remote Microsoft desktop to connect.     

2. For the code to run in its entirety, it will take hours to days possibly (based on the size of the video data being processed and how the data is being processed). To avoid running into small mistakes that will force you to run the code all over again, execute this program section by section and verify the processing worked correctly by looking at the supplemental figures. 

3. This code has paths to directories and files written using the Window's forward slash convention. This may impact Mac users using this code. 

## Additional notes: 

1. Vignetting is defined in photography as the reduction of an image's brightness or saturation toward the periphery compared to the image center. For the video camera, the sun glint off the ocean surface causes natural vignetting where the image is brightest in the center of the sun glint and saturated elsewhere. This makes it challenging to identify wave breaking in the images.
2. Images are georeferenced using the Trimble software where each pixel in the image is assigned a lat/lon or UTM coordinate based on the position/attitude of the aircraft and the distortion of the lens. Additionally, boresight adjustments (i.e., offsets between reference frames the GPS-IMU, and the video camera) are accounted for.
