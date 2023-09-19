# Main.m Script

Script for estimating lambda of c and other wave-breaking statistics 
from processed video data using the IO Industries FLARE 12M125-CL
camera (4096 by 3072-pixel resolution, 10-bit, 5-Hz sampling rate) 
collected onboard a research aircraft. Images are georeferenced using a 
Novatel SPAN LN200 Inertial Motion Unit (IMU) coupled to a ProPak6 GPS
receiver which provided aircraft trajectory information 
(attitude and position).

This script imports the following processed data: 

1. Video camera images georeferenced with Trimble software: Consists of GeoTIFF files where each pixel in the image is assigned a lat/lon or UTM coordinate based on the position/attitude of the aircraft and the distortion of the lens. Additionally, boresight adjustments (i.e., offsets between reference frames the GPS-IMU, and the video camera) are accounted for. 

2. Aircraft trajectory and attitude: Consists of the EO text files that provide the position and attitude of the aircraft over time from the coupled GPS-IMU onboard the research aircraft. This time series is interpolated to match the time steps of the video camera images.  

For specifics on video and trajectory data processing, see Nick Statom's
documentation [here](https://docs.google.com/document/d/1qbaBH98IW1tJrMfxC6TKL-jQIcMPm7KK_KRrxk3rBb0/edit).

The lambda of c calculation steps include the following: 

1. Load aircraft 

Below are a few tips for running the code. 

1. This code can be run locally on your computer or remotely through whichever airseaserver the data is located on. Either way, your local computer needs to be connected to the UCSD VPN; use Cisco AnyConnect to do this. If running the code locally, connect to the remote server  through your file system explorer to gain remote access to the data (i.e., mount to the remote server on your local computer). If running the code on the airseaserver, use a remote Microsoft desktop to connect.     

2. For the code to run in its entirety, it will take hours to days possibly (based on the size of the video data being processed and how the data is being processed). To avoid running into small mistakes that will force you to run the code all over again, execute this program section by section and verify the processing worked correctly by looking at the supplemental figures. 

3. This code has paths to directories and files written using the Window's forward slash convention. This may impact Mac users using this code. 

Some additional notes: 

1. Vignetting is defined in photography as the reduction of an image's brightness or saturation toward the periphery compared to the image center. For the video camera, the sun glint off the ocean surface causes natural vignetting where the image is brightest in the center of the sun glint and saturated elsewhere. This makes it challenging to identify wave breaking in the images.
