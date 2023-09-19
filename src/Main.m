%% Main Script for estimating Lambda of c distributions. 
% Authors: Teodor Vrecica (tvrecica@ucsd.edu) and Luke Colosi (lcolosi@ucsd.edu)

%--------------------------------- Notes ---------------------------------%
% Script for estimating lambda of c and other wave breaking statistics 
% from processed video data using the IO Industries FLARE 12M125-CL
% camera (4096 by 3072 pixel resolution, 10 bit, 5-Hz sampling rate) 
% collected onboard a research aircraft. Images are georeferenced using a 
% Novatel SPAN LN200 Inertial Motion Unit (IMU) coupled to a ProPak6 GPS
% reciever which provided aircraft trajectory information 
% (attitude and postion).
% 
% This script imports the following processed data: 
% 
% (1) Video camera images georeferenced with Trimble software: Consist of
%     geotiff files where each pixel in the image is assigned a lat/lon or
%     UTM coordinate based on the position/atitude of the aircraft and the
%     distortion ofthe lens. Additionally, boresight adjustments (i.e.,
%     offsets between reference frames the GPS-IMU and the video camera)
%     are accounted for. 
%
% (2) Aircraft trajectory and attitude: Consist of the EO text files that
%     provide the position and attitude of the aircraft over time from the 
%     coupled GPS-IMU onboard the research aircraft. This time series is
%     interpolated to match the time steps of the video camera images.  
% 
% For specifics on video and trajectory data processing, see Nick Statom's
% documentation here: 
% 
% https://docs.google.com/document/d/1qbaBH98IW1tJrMfxC6TKL-jQIcMPm7KK_KRrxk3rBb0/edit
% 
% The lambda of c calculation steps include the following: 
% 
% (1) Load aircraft 
% 
% Below are a few tips for running the code. 
% 
% (1) This code can be ran locally on your computer or remotely through
%     whichever airseaserver the data is located on. Either way, your local
%     computer needs to be connected to the ucsd vpn; use cisco AnyConnect
%     to do this. If running the code locally, connect to remote server 
%     through your file system explorer to gain remote access to the data 
%     (i.e., mount to the remote server on your local computer). 
%     If running the code on the airseaserver, use remote microsoft desktop
%     to connect.     
% 
% (2) For the code to run in its entirety, it will take hours to days 
%     possibly (based on the size of the video data being processed and 
%     how the data is being processed). To avoid running into a small 
%     mistake which will force you to run the code all over again, execute
%     this program section by section and verify the processing worked
%     correctly by looking at the supplemental figures. 
% 
% (3) This code has paths to directories and files written using the
%     windows forward slash convention. This may impact mac users using
%     this code. 
% 
% Some additional notes: 
% 
% (1) Vignetting is definited in photography as the reduction of an 
%     image's brightness or saturation toward the periphery compared to the
%     image center. For the video camera, the sun glint off the ocean 
%     surface causes natural vignetting where the image is brightest in
%     the center of the sun glit and saturated elsewhere. This makes it
%     challenging to identify wave breaking in the images.
% 
%-------------------------------------------------------------------------%

clc, clearvars, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- This is the only section of code you will need to change --------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set path to raw (non-georeferenced) video images for a given flight
dirRaw = 'X:\\TFO_2021\Processed\VIDEO\16bit_TIF_Frames\20210519\';

% Set path to georeferenced video images for a given flight 
dirout = 'X:\\TFO_2021\Processed\VIDEO\Trimble\20210519\';

% Set path to save verification figures
dirV = 'D:\DEPLOYMENTS\TFO_2021\figs\video\Quality_Control\';

% Set flight stability criteria parameters
maxPer=[25 25];                                                             % Maximum percent of flight track to be removed at the begin and end of the track. For example, maxPer = [25,25] means that at a minimium, 50 percent(from 25% - 75%) of the flight track will be used in analysis with the 0-25% and 75%-100% removed from the flight track.  
sigRoll=3;                                                                  % Maximum allowed roll standard deviations for a stable segment.    
sigPitch=3;                                                                 % Maximum allowed pitch standard deviations for a stable segment.
sigHeading=6;                                                               % Maximum allowed heading standard deviations for a stable segment.

% Specify whether verification figures will be plotted
ploton = true; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a directory for verification figures if one does not exist
if ~exist(dirV, 'dir'), mkdir(dirV); end

%% Load EO data and determine time indices of each track in the EO file

% Rename EO directory in dirout path as Trimble_EO 
if isfolder([dirout 'EO'])
    movefile([dirout 'EO'],[dirout 'Trimble_EO']);
end

% Obtain filenames of all EO files 
D=dir([dirout 'Trimble_EO\' '*EO.txt']);

% Set path to EO file and its directory  
EOpath=[dirout 'Trimble_EO\' D(1).name];
EOdir=[dirout 'Trimble_EO\'];

% Obtain number of lines in the EO text file header 
Counter=HeaderPosition(EOpath,'(sec)');

% Load trajectory and attitude data 
A = importdata(EOpath,' ',Counter);

% Get the text component of data (timestamp into column 1, track number and image number into column 2)
An=cropText(A);

% Determine the start and end time indices for each flight tracks.
tracks = DefineTracks(An);

%% Determine aircraft stability for each flight track 
% Note: If flight track is NOT steady, Lambda of c distributions and 
% whitecap coverage will be determined.  

% Identify track stability based on input criteria
trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,An);

% Plot out headings in flight tracks, with all other data (used portion of 
% the track, whether is it labeled as steady or not).
if ploton == true

    % Create a subdirectory for flight trajectory and stability plots
    dirstab=[dirV 'Flight_Trajectories_And_Stability\'];
    if ~exist(dirstab, 'dir'), mkdir(dirstab); end
    
    % Generate plots
    plotOutTracks(A,tracks,trackTag,dirstab,An,dirout);
end

%% Remove vignetting and crop out high glint 

% Obtain the filenames of the non-georeferenced video images 
D_Im=dir([dirRaw '*.tif']);

% Define tracks in images
tracks_Im = DefineTracksIm(D_Im,tracks);

return
% First we obtain the mean of the image and of std for each track. This is
% computed for the stable portion of the tracks, as indicated in trackTag
% variable.

% The RM_Nr and meanOriginal variables check for significant variations in
% image lighting. RM_Nr represents image numbers to be removed, and
% meanOriginal is the mean lighting in every image

winSize=7;
[meanIm,stdIm,RM_Nr,meanOriginal]=getImMeanStd_V5(dirRaw,D_Im,tracks_Im,trackTag,tracks,winSize,dirout);

dirV=[dirout 'Saved_Matlab_Data\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end
save([dirV 'mean_Images_Per_Track.mat'],'meanIm','stdIm','RM_Nr','meanOriginal');


% Next, the raw images are processed (std and mean brightness are
% equalized), and are saved to the Output folder (dirout\Output).


% Export images with vignetting removed
% dirProcessed='\\Airseaserver28\d\MASS\RAW\LCDRI_2017\VIDEO\TIFF_16BIT_EXPORT_Corrected\20170323\';
dirProcessed=[dirout 'nongeoreferenced_Images_With_Vignetting_Removed_V2\'];
if ~exist(dirProcessed, 'dir'), mkdir(dirProcessed); end
mina2=saveImage_test(dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirProcessed);

%Save mask images
dirProcessed=[dirout 'mask_For_Determining_High_Glint_Areas_V2\'];
if ~exist(dirProcessed, 'dir'), mkdir(dirProcessed); end
mina2=saveImageMask(dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirProcessed);

% Eliminate areas with high glint. Std mag determines the number of
% standard deviations for categorizing high glint areas.
stdMag=-0.1;
[Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag);

dirV=[dirout 'Quality_Control\Vignetting\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end
% Finally, the original and processed images are plotted out for
% comparison. The mean image and mean std are also plotted out. 
plotProcessedImages(dirV,dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,Glint);


%% Build and execute .bat file for georeferencing processed images


% Split the EO files into individual tracks (think file organization is
% better this way).
cleanEOfile(EOpath,EOdir,tracks,D);

% A matlab script which builds .bat file for georeferencing images is
% called. 
camName = 'flare'; %current camera options are 'jai', 'flare', 'sc6000', 'sc6700', 'd810', 'doppvis'
utmZone='10 N';
%% !!!!!!!!!!!!!!!!!!!!!!!!!! added 10 to resolution
for i=11:12%1:length(tracks)
    if trackTag(i).stable==1
        xx=nanmean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1));
        yy=nanmean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2));
        [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
        
        dirin_img = [dirout 'nongeoreferenced_Images_With_Vignetting_Removed_V2\Track_' num2str(i) '\']; %directory of all tif imagery to use for this single flight
        dirin_prj = dirout;
        
        buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,i,camName);
    end
end

% Build .bat file for the mask
for i=11:12%1:length(tracks)
    if trackTag(i).stable==1
        xx=nanmean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1));
        yy=nanmean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2));
        [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
        
        dirin_img = [dirout 'mask_For_Determining_High_Glint_Areas\Track_' num2str(i) '\']; %directory of all tif imagery to use for this single flight
        dirin_prj = dirout;
        
        buildBatMask(dirin_img,dirin_prj,proj_lat,proj_lon,i,camName);
    end
end

% Execute .bat files
for i=11:12%1:length(tracks)
    if trackTag(i).stable==1
        fileNN=dir([dirout 'Project\Track_' num2str(i) '\' '*.bat']);
        system([dirout 'Project\Track_' num2str(i) '\' fileNN.name]);
%         while ~exist([dirout 'Output\Track_' num2str(i) '\blank.txt'])
%            pause(300);
%            fprintf('ne')
%         end
        pause(150);
    end
end

for i=11:12%1:length(tracks)
    if trackTag(i).stable==1
        fileNN=dir([dirout 'Project\TrackMask_' num2str(i) '\' '*.bat']);
        system([dirout 'Project\TrackMask_' num2str(i) '\' fileNN.name]);
%         while ~exist([dirout 'Output\TrackMask_' num2str(i) '\blank.txt'])
%            pause(300);
%            fprintf('ne')
%         end
        pause(150);
    end
end


%% Determine brightness threshold for each track
GlobalOrLocal=0; % Check if a global brightness threshold is defined for the area (1), or a local one for each image (0).
arrayHist=0:16:256*256;% For 16bit
dirProcessed=[dirout 'Output\'];
localStep=3;%A range determining how many images are included into the analysis for local threshold (200 would indicate up to 100 images taken before or after the current one)
peakPercentage=20;%what percentage of peak magnitude is taken as a breaking threshold (i.e. 10 would mean 10%)

if GlobalOrLocal==1
    [Threshold,N]=global_Thresh(peakPercentage,dirProcessed,D_Im,tracks_Im,trackTag,tracks,arrayHist,meanIm,RM_Nr,Glint,Glint_mask);
elseif GlobalOrLocal==0
    [Threshold,N]=local_Thresh(peakPercentage,dirProcessed,D_Im,tracks_Im,trackTag,tracks,arrayHist,meanIm,RM_Nr,Glint,Glint_mask,localStep);
end

for i=11:12
    N2=mean(N(i).hi);
    N_CS=cumsum(N2(200:500));
    N_CS=N_CS/nanmax(N_CS);
    N_CS=1-N_CS;
    ind1=find(N_CS<0.5,1,'first');
    [res_x, idx_of_result]=knee_pt(N_CS(ind1:301),ind1:301);
    Threshold(i).th(:)=10000%(199+ind1+idx_of_result)*16;
end

for i=11:12
    for j=1:length(N(i).hi(:,1))
        if N(i).hi(j,1)==0
            Threshold(i).th(j)=0;
        else
            N2=(N(i).hi(j,:));
            N_CS=cumsum(N2(200:500));
            N_CS=N_CS/nanmax(N_CS);
            N_CS=1-N_CS;
            ind1=find(N_CS<0.5,1,'first');
            [res_x, idx_of_result]=knee_pt(N_CS(ind1:301),ind1:301);
            Threshold(i).th(j)=(199+ind1+idx_of_result)*16;
        end
        
    end
end

if GlobalOrLocal==1
    dirV=[dirout 'Saved_Matlab_Data\'];
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    save([dirV 'thresh_global.mat'],'Threshold','N');
elseif GlobalOrLocal==0
    Threshold=denoiseThreshold(Threshold,tracks,trackTag);
    dirV=[dirout 'Saved_Matlab_Data\'];
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    save([dirV 'thresh_local_NewMethod.mat'],'Threshold','N');
end
% plot out graphs for determining threshold, and give some examples
dirV=[dirout 'Quality_Control\Brightness_Threshold_New\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end
plot_Thresh(Threshold,N,dirProcessed,dirV,D_Im,tracks_Im,trackTag,tracks,arrayHist,Glint,Glint_mask,meanIm,GlobalOrLocal);

%% temp code for thresholding


%% Start running the code for determining Lambda and wc coverage
%new - 0 
%new2 - 10
%new3 - 5
%new4 - 15
%

% v3 - portion of image v4- portion, 30cm v5 - old processing portion of
% image
for i=11:1:12%1:length(tracks)
    if trackTag(i).stable==1
        if i>=10
            dirV=[dirout 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_' num2str(i) '\'];
        else
            dirV=[dirout 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_0' num2str(i) '\'];
        end
        if ~exist(dirV, 'dir'), mkdir(dirV); end
        
       
        
        dirinIM=[dirout 'Output\Track_' num2str(i) '\'];
        dirinIM_mask=[dirout 'Output\TrackMask_' num2str(i) '\'];
        filenamesIM=dir([dirinIM '*.tif']);
        
        offSetEnds=0;
        
        i1=1+3+(trackTag(i).range(1)+offSetEnds)-tracks(i).Indices(1);
        i2=1+(trackTag(i).range(2)-offSetEnds)-tracks(i).Indices(1)-1;
        %Threshold(i)=11000;
        %process_Lambda_New(i1,i2,filenamesIM,dirinIM,dirV,Threshold(i),300,0.1,0.2)
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(2);
        strN=filenamesIM(1).name;
        dxIm = str2num(strN(pozicija+1:pozicija+2))/100;
%         dxIm = 0.3;
        % Define angles (from mean wind direction) which are considered as
        % breaking.
        SpreadMax=110;
        
        % Determine wind direction from the data. The images are indexed in
        % seconds from the begining of the week. This requires entering the
        % start date manually.
        StartTime=datenum([2021,5,2,0,0,0]);
        % Note - some days lack data in SONIC2, but it is prefered option
        % otherwise
        load('PLD2_1hz_ALL_1.mat');
        WindTime=PLD2_1hz.time;
        Speeds=PLD2_1hz.wind_direction;
        Speeds=PLD2_1hz.TWD;
        WindDir=determine_wind(Speeds,WindTime,StartTime,tracks(i).Indices,An);
        process_Lambda_New(i1,i2-3,filenamesIM,dirinIM,dirV,Threshold(i).th,300,dxIm,0.2,WindDir,SpreadMax,Glint(i).list,dirinIM_mask,meanIm(i).im,RM_Nr(i).nr,GlobalOrLocal,trackTag(i),tracks(i))
    end
end



% Plot out randomly chosen examples from the track
dirV=[dirout 'Quality_Control\LambdaIms_inbuilt\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end

dirV2=[dirout 'Plots\DistributionsAll\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirout 'Plots\DistributionsAll_WC_Glint\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirout 'Plots\DistributionsAll_new\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirout 'Plots\Distributions_Inbuilt\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end


dirV3=[dirout 'Plots\Summary\'];
if ~exist(dirV3, 'dir'), mkdir(dirV3); end

dirV4=[dirout 'Plots\Alongtrack_OppositeDirection\'];
if ~exist(dirV4, 'dir'), mkdir(dirV4); end

dirV4=[dirout 'Plots\Alongtrack_test3\'];
if ~exist(dirV4, 'dir'), mkdir(dirV4); end
addpath('utilities')


dirV2=[dirout 'Plots\DistributionsAll_new\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end
%new - 0 
%new2 - 10
%new3 - 5
%new4 - 15
for i=11:1:11%1:length(tracks)
    i
    if trackTag(i).stable==1
        if i>=10
            dirL=[dirout 'Lambdas_inbuilt\Track_' num2str(i) '\'];
        else
            dirL=[dirout 'Lambdas_inbuilt\Track_0' num2str(i) '\'];
        end
        filenames=dir([dirL '*.mat']);
        dirinIM=[dirout 'Output\Track_' num2str(i) '\'];
        
        filenamesIM=dir([dirinIM '*.tif']);
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(2);
        strN=filenamesIM(1).name;
        dxIm = str2num(strN(pozicija+1:pozicija+2))/100;
        
        PARAM.dt = 0.2;
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(3);
        strN=filenamesIM(1).name;
        PARAM.dx = str2num(strN(pozicija+1:pozicija+2))/100;
        
        dirin_img = [dirout 'Output\Track_' num2str(i) '\'];
        filenamesI=dir([dirin_img '*.tif']);
        
        %determineDists(dirV2,i,dirL,filenames,dxIm,tracks,A,Threshold);
        i1=1+3+trackTag(i).range(1)-tracks(i).Indices(1);
        i2=1+trackTag(i).range(2)-tracks(i).Indices(1)-1-3;
        %plotSummary(dirV3,i,dirL,filenames,dxIm,tracks,A,dirProcessed,i1,i2,RM_Nr(i).nr);
        
        PlotExamples(dirV,i,tracks,dirL,filenames,dirin_img,filenamesI,PARAM);
        
        %determineDists_Along(dirV4,i,dirL,filenames,dxIm,tracks,A,dirinIM);
    end
end
windData='PLD2_1hz.mat';
waveData='WG1mc_Snn_Hs_30min.mat';

% determine params for tfo
for i=2:3%length(trackTag)
    if trackTag(i).stable==1
        
        
        WGTime=PLD2_AVG.time;
        WGTime2=PLD2_FFT.time;
        Hs=PLD2_AVG.Hs_fft;
        Freq=0:0.01:0.82;
        Spectra=[];
        Spectra2=[];
        for ii=1:83
            Spectra(ii,:)=eval(['PLD2_FFT.avgSzz_' num2str(ii)]);
            Spectra2(ii,:)=eval(['PLD2_FFT.avgSzz_' num2str(ii)]);
        end
        [Hs_all(i),freq_all(i)]=determine_params_TFO(Spectra,Freq,Hs,WGTime,WGTime2,StartTime,tracks(i).Indices,An,Spectra2);
        
        WindTime=PLD2_1hz.time;
        us=movmean(PLD2_1hz.TWS*0.41/log(1.11/0.0002),100);
        ustar(i)=determine_ustar(us,WindTime,StartTime,tracks(i).Indices,An);
    end
end

plotAll_TFO(dirV2,trackTag,StartTime,tracks,An,Hs_all,freq_all,ustar);

% Plot out other verification plots (wc coverage in images, total length of
% breaking, spectrogram of Lambda distribution).

%plotDists();

%% Development Code

% Load buoy positions
% load('D:\MASS\Processed\L_Computations\L_Computations\Test\TFOEx21_DEP01_BuoyGPS.mat')











