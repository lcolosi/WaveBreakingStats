%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Script for estimating Lambda of c distributions  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See Github repo (https://github.com/lcolosi/WaveBreakingStats/src/README.md) 
% for documentation.

clc, clearvars, close all

display_text('Estimating Lambda of c distributions.','title')
display_text('Step 1: Setting input parameters.','section')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- This is the only section of code you will need to change --------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select process to run (0 or 1)
option_plot          = 0; 
option_eo_split      = 0;
option_image_proc    = 0;                                                   
option_globalOrlocal = 1;                                                 

% Start date of flight in UTC time
StartDate = '20210519';                                                    

% Directories
dirRaw = 'X:\\TFO_2021\Processed\VIDEO\16bit_TIF_Frames\20210519\';         
dirProc = 'X:\\TFO_2021\Processed\VIDEO\Trimble\20210519\';                 
dirV = 'D:\DEPLOYMENTS\TFO_2021\figs\video\Quality_Control\20210519\';      
dirOut = 'D:\DEPLOYMENTS\TFO_2021\data\video\Intermediate_Products\20210519\'; 

% Flight stability criteria parameters
maxPer = [25 25];                                                           
sigRoll = 3;                                                                   
sigPitch = 3;                                                               
sigHeading = 6;                                                             
Shift = [-0.8*maxPer(1) 0.8*maxPer(2);0.8*maxPer(1) -0.8*maxPer(2);0 0];    
Nstd = 4;                                                                   
tCheck = 7;                                                                 

% Vignette removal parameters
winSize = 7;                                                                
sigma_ff_m = 300;                                                           
sigma_ff_v = 450;                                                           
sigBrightness = 3;                                                          
B_threshold = 0.8;                                                          
n_sigma = 5;                                                                
stdMag=-0.1;                                                                

% Trimble georeferencing project parameters 
camName = 'flare';                                                         
utmZone='10 N';                                                            

% Brightness threshold parameters
localStep=3;                                                               
peakPercentage=20;                                                         

% Set text interpreter
set_interpreter('latex')

display_text(['Date: ' StartDate(end-3:end-2) '/' StartDate(end-1:end) '/' StartDate(1:4)],'body')
display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import EO data and determine time indices of each track in the EO file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

display_text('Step 2: Loading EO data and identifing flight tracks.','section')

% Identify that EO files are processed through Trimble
if isfolder([dirProc 'EO'])
    movefile([dirProc 'EO'],[dirProc 'Trimble_EO']);
end

% Obtain filenames of all EO files 
D = dir([dirProc 'Trimble_EO\' '*EO.txt']);

% Set path to EO file and its directory  
EOpath = [dirProc 'Trimble_EO\' D(1).name];
EOdir = [dirProc 'Trimble_EO\'];

% Obtain number of lines in the EO text file header 
Counter = HeaderPosition(EOpath,'(sec)');

% Load trajectory and attitude data 
A = importdata(EOpath,' ',Counter);

% Get the text component of data (timestamp into column 1, track number and image number into column 2)
An = cropText(A);

% Determine the start and end time indices for each flight tracks.
tracks = DefineTracks(An);

% Split the EO files into individual track files (helps to break up code
% when georeferencing images with Trimble)
if option_eo_split == 1
    cleanEOfile(EOpath,EOdir,tracks,D);
end

% Convert from gps time in seconds of the week to UTC time
gps_time = str2num(cell2mat(An(:,1)));                                      %#ok<ST2NM> 
utc_time = gpssw2utcdn(gps_time,StartDate);

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine aircraft stability for each flight track 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Note: If flight track is NOT steady, Lambda of c distributions and 
% whitecap coverage will not be computed.  

display_text('Step 3: Determining periods of aircraft stability for each track.','section')

% Identify track stability based on input criteria
trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,Shift,Nstd,tCheck,An);

% Plot headings in flight tracks, with all other data (used portion of 
% the track, whether is it labeled as steady or not).
if option_plot == 1

    % Create a subdirectory for flight trajectory and stability plots
    dirTS = [dirV 'Flight_Trajectories_And_Stability\'];
    if ~exist(dirTS, 'dir'), mkdir(dirTS); end
    
    % Set up stability criteria matrix
    sc = [sigRoll,sigPitch,sigHeading,Nstd];

    % Generate plots
    plotOutTracks(A,tracks,trackTag,dirTS,An,dirProc,sc,utc_time,StartDate);
end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove vignetting and computing parameters for identifying high glint regions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_text('Step 4: Removing vignetting, equalizing image, and computing images for identifying high sun glint regions.','section')

% Obtain the filenames of the non-georeferenced video images 
D_Im = dir([dirRaw '*.tif']);

% Define tracks in images
tracks_Im = DefineTracksIm(D_Im,tracks);

% Preform image processing of raw non-georeference images
if option_image_proc == 1

    % Create a subdirectory for intermediate image products
    dirVR=[dirOut 'nongeoreferenced_images_with_Vignetting_Removed\'];
    if ~exist(dirVR, 'dir'), mkdir(dirVR); end
    dirMGlint=[dirOut 'mask_for_Determining_High_Glint_Areas\'];
    if ~exist(dirMGlint, 'dir'), mkdir(dirMGlint); end

    % Compute the mean brightness and standard deviation of each track
    [meanIm,stdIm,RM_Nr,meanOriginal]=getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,tracks,winSize,sigma_ff_m,B_threshold,n_sigma);
    
    % Save output from getImMeanStd 
    save([dirOut 'mean_Images_Per_Track.mat'],'meanIm','stdIm','RM_Nr','meanOriginal');
    
    % Remove vignetting and equalize raw nongeoreferenced image; save 
    % intermediate products to dirVR
    saveImage(dirRaw,D_Im,tracks_Im,trackTag,tracks,dirVR,sigma_ff_v,n_sigma,meanOriginal,sigBrightness);
    
    % Save normalized image for detecting high sun glint areas; save 
    % intermediate products to dirMGlint
    saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,dirMGlint);
    
    % Eliminate areas with high glint. Std mag determines the number of
    % standard deviations for categorizing high glint areas.
    [Glint]=det_glint(meanIm,stdMag,tracks,trackTag);

    % Save output from det_glint 
    save([dirOut 'glint_threshold_Per_Track.mat'],'Glint','Glint_mask');

else 
    
    % Load mean brightness and standard deviation imagesand glint threshold
    load([dirOut 'mean_Images_Per_Track.mat'])
    load([dirOut 'glint_threshold_Per_Track.mat'])

end

% Plot 
if option_plot == 1

    % Create a subdirectory for vignette removal plots
    dirVn=[dirV 'Vignetting\'];
    if ~exist(dirVn, 'dir'), mkdir(dirVn); end

    % Generate plots 
    plotProcessedImages(dirVn,dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,Glint);
end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build and execute .bat file for georeferencing processed images with Trimble 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

display_text('Step 5: Building and executing .bat files for georeferencing processed images','section')

%--- Build batch file for georeferencing images ---%

% Loop through tracks
for i=1:length(tracks)

    % Check if track is stable
    if trackTag(i).stable==1

        % Compute mean northing and easting UTM coordinate of ith track and
        % convert it to decimal degree longitude and latitude
        xx=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1),'omitnan');
        yy=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2),'omitnan');
        [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
        
        % Set the paths to nongeoreferenced images and the output directory
        % for post-processed georeferenced images 
        dirin_img = [dirProc 'nongeoreferenced_Images_With_Vignetting_Removed\Track_' num2str(i) '\']; 
        dirin_prj = dirProc;
        
        % Build batch file for running trimble georeferencing project 
        buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,i,camName);
    end
end

%--- Build batch file for georeferencing images for identifying high glint regions ---%

% Loop through tracks 
for i=1:length(tracks)

    % Check if track is stable
    if trackTag(i).stable==1
        
        % Compute mean northing and easting UTM coordinate of ith track and
        % convert it to decimal degree longitude and latitude
        xx=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1),'omitnan');
        yy=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2),'omitnan');
        [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
        
        % Set the paths to masks for determing high glint and output
        % directory for georeferenced masks
        dirin_img = [dirProc 'mask_For_Determining_High_Glint_Areas\Track_' num2str(i) '\'];
        dirin_prj = dirProc;
        
        % Build batch file for running trimble georeferencing project
        buildBatMask(dirin_img,dirin_prj,proj_lat,proj_lon,i,camName);
    end
end

%--- Execute batch files ---%

% Loop through tracks 
for i=1:length(tracks)
    
    % Check if track is stable 
    if trackTag(i).stable==1
        
        % Obtain file name of .bat file  
        fileNN=dir([dirProc 'Project\Track_' num2str(i) '\' '*.bat']);
        
        % Execute batch file for trimble georeferencing image project
        system([dirProc 'Project\Track_' num2str(i) '\' fileNN.name]);

    end
end

% Loop through tracks
for i=1:length(tracks)
    
    % Check if track is stable 
    if trackTag(i).stable==1
        
        % Obtain file name of .bat file
        fileNN=dir([dirProc 'Project\TrackMask_' num2str(i) '\' '*.bat']);
        
        % Execute batch file for trimble georeferencing mask project
        system([dirProc 'Project\TrackMask_' num2str(i) '\' fileNN.name]);

    end
end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine brightness threshold of breakers for each track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set pixel value range for 16 bit images
arrayHist=0:16:256*256;

% Set directory to save brightness threshold
dirProcessed=[dirProc 'Output\'];

% Jessica's Method
if option_globalOrlocal == 1
    [Threshold,N]=global_Thresh(peakPercentage,dirProcessed,D_Im,tracks_Im,trackTag,tracks,arrayHist,meanIm,RM_Nr,Glint,Glint_mask);
else 
    [Threshold,N]=local_Thresh(peakPercentage,dirProcessed,D_Im,tracks_Im,trackTag,tracks,arrayHist,meanIm,RM_Nr,Glint,Glint_mask,localStep);
end

% Slight variation on Jessica's method

% Loop through tracks
for i=1:length(tracks)
    for j=1:length(N(i).hi(:,1))
        if N(i).hi(j,1)==0
            Threshold(i).th(j)=0;
        else
            N2=(N(i).hi(j,:));
            N_CS=cumsum(N2(200:500));
            N_CS=N_CS/max(N_CS,[],'omitnan');
            N_CS=1-N_CS;
            ind1=find(N_CS<0.5,1,'first');
            [res_x, idx_of_result]=knee_pt(N_CS(ind1:301),ind1:301);
            Threshold(i).th(j)=(199+ind1+idx_of_result)*16;
        end
        
    end
end

% 
if option_globalOrlocal == 1
    dirV=[dirProc 'Saved_Matlab_Data\'];
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    save([dirV 'thresh_global.mat'],'Threshold','N');
else 
    Threshold=denoiseThreshold(Threshold,tracks,trackTag);
    dirV=[dirProc 'Saved_Matlab_Data\'];
    if ~exist(dirV, 'dir'), mkdir(dirV); end
    save([dirV 'thresh_local_NewMethod.mat'],'Threshold','N');
end

% Plot 
if option_plot

    % Create a subdirectory for vignette removal plots
    dirBT=[dirV 'Quality_Control\Brightness_Threshold_New\'];
    if ~exist(dirBT, 'dir'), mkdir(dirBT); end

    % Generate plots
    plot_Thresh(Threshold,N,dirProcessed,dirV,D_Im,tracks_Im,trackTag,tracks,arrayHist,Glint,Glint_mask,meanIm,GlobalOrLocal);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start running the code for determining Lambda and white cap coverage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v3 - portion of image v4- portion, 30cm v5 - old processing portion of
% image
for i=11:1:12%1:length(tracks)
    if trackTag(i).stable==1
        if i>=10
            dirV=[dirProc 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_' num2str(i) '\'];
        else
            dirV=[dirProc 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_0' num2str(i) '\'];
        end
        if ~exist(dirV, 'dir'), mkdir(dirV); end
        
       
        
        dirinIM=[dirProc 'Output\Track_' num2str(i) '\'];
        dirinIM_mask=[dirProc 'Output\TrackMask_' num2str(i) '\'];
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
dirV=[dirProc 'Quality_Control\LambdaIms_inbuilt\'];
if ~exist(dirV, 'dir'), mkdir(dirV); end

dirV2=[dirProc 'Plots\DistributionsAll\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirProc 'Plots\DistributionsAll_WC_Glint\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirProc 'Plots\DistributionsAll_new\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end

dirV2=[dirProc 'Plots\Distributions_Inbuilt\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end


dirV3=[dirProc 'Plots\Summary\'];
if ~exist(dirV3, 'dir'), mkdir(dirV3); end

dirV4=[dirProc 'Plots\Alongtrack_OppositeDirection\'];
if ~exist(dirV4, 'dir'), mkdir(dirV4); end

dirV4=[dirProc 'Plots\Alongtrack_test3\'];
if ~exist(dirV4, 'dir'), mkdir(dirV4); end
addpath('utilities')


dirV2=[dirProc 'Plots\DistributionsAll_new\'];
if ~exist(dirV2, 'dir'), mkdir(dirV2); end
%new - 0 
%new2 - 10
%new3 - 5
%new4 - 15
for i=11:1:11%1:length(tracks)
    i
    if trackTag(i).stable==1
        if i>=10
            dirL=[dirProc 'Lambdas_inbuilt\Track_' num2str(i) '\'];
        else
            dirL=[dirProc 'Lambdas_inbuilt\Track_0' num2str(i) '\'];
        end
        filenames=dir([dirL '*.mat']);
        dirinIM=[dirProc 'Output\Track_' num2str(i) '\'];
        
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
        
        dirin_img = [dirProc 'Output\Track_' num2str(i) '\'];
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











