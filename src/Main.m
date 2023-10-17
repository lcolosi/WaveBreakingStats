%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Script for estimating Lambda of c distributions  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See Github repo (https://github.com/lcolosi/WaveBreakingStats/src/README.md) 
% for documentation.

clc, clearvars -except meanIm stdIm RM_Nr meanOriginal Glint , close all

display_text('Estimating Lambda of c distributions.','title')
display_text('Step 1: Setting input parameters.','section')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- This is the only section of code you will need to change --------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select process to run (0 or 1)
option_plot          = 1; 
option_eo_split      = 1;
option_image_proc    = 1;  
option_run_bat       = 1;
option_globalOrlocal = 1;                                                 

% Experiment title
exp = 'SMODE IOP1';

% Start date of flight in UTC time
StartDate = '20221025';                                                    

% Directories
dirRaw = 'W:\\SMODE_2022\RAW\VIDEO\Frames\20221025_1\Images\';         
dirProc = 'W:\\SMODE_2022\RAW\VIDEO\Frames\20221025_1\';                  
dirV = 'W:\\SMODE_2022\RAW\VIDEO\Frames\20221025_1\QC_Plots\';      
dirOut = 'W:\\SMODE_2022\PROCESSED\VIDEO\20221025_1\';

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
sigBrightness = 3000;                                                          
B_threshold = 0.8;                                                          
n_sigma = 5;                                                                
stdMag=-0.1;                                                                

% Trimble georeferencing project parameters  
res = 5;                        
utmZone='10 N';                                                            

% Brightness threshold parameters
localStep=3;                                                               
peakPercentage=20;

% Lambda of c parameters
SpreadMax=110;

% Set text interpreter
set_interpreter('latex')

display_text(['Experiment: ' exp],'body')
display_text(['Date: ' StartDate(end-3:end-2) '/' StartDate(end-1:end) '/' StartDate(1:4)],'body')
display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import EO data and determine time indices of each track in the EO file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

display_text('Step 2: Importing EO data and identifying flight tracks.','section')

% Obtain filenames of all EO files 
D = dir([dirProc 'EO\' '*EO.txt']);

% Set path to EO file and its directory  
EOpath = [dirProc 'EO\' D(1).name];
EOdir = [dirProc 'EO\'];

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

display_text('Step 3: Determining periods of aircraft stability.','section')

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
    plotOutTracks(A,tracks,trackTag,dirTS,An,sc,utc_time,StartDate,utmZone);
end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove vignetting and computing parameters for identifying high glint regions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_text('Step 4: Removing vignetting, equalizing image, and generating images for identifying high sun glint regions.','section')

% Obtain the filenames of the non-georeferenced video images 
D_Im = dir([dirRaw '*.tif']);

% Define tracks in images
tracks_Im = DefineTracksIm(D_Im,tracks);

% Preform image processing of raw non-georeference images
if option_image_proc == 1

    % Create a subdirectory for intermediate image products
    dirVR=[dirProc '\Images\Intermediate_Products\nongeoreferenced_images_with_Vignetting_Removed\'];
    if ~exist(dirVR, 'dir'), mkdir(dirVR); end
    dirMGlint=[dirProc '\Images\Intermediate_Products\mask_for_Determining_High_Glint_Areas\'];
    if ~exist(dirMGlint, 'dir'), mkdir(dirMGlint); end

    % Compute the mean brightness and standard deviation of each track
    [meanIm,stdIm,RM_Nr,meanOriginal]=getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,tracks,winSize,sigma_ff_m,B_threshold,n_sigma);
    
    % Save output from getImMeanStd 
    save([dirProc '\Images\Intermediate_Products\mean_Images_Per_Track.mat'],'meanIm','stdIm','RM_Nr','meanOriginal');
    
    % Remove vignetting and equalize raw nongeoreferenced image; save 
    % intermediate products to dirVR
    saveImage(dirRaw,D_Im,tracks_Im,trackTag,tracks,dirVR,sigma_ff_v,n_sigma,meanIm,stdIm,meanOriginal,sigBrightness);
    
    % Save normalized image for detecting high sun glint areas; save 
    % intermediate products to dirMGlint
    saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,dirMGlint);
    
    % Eliminate areas with high glint. StdMag determines the number of
    % standard deviations for categorizing high glint areas.
    [Glint]=det_glint(meanIm,stdMag,tracks,trackTag);

    % Save output from det_glint 
    save([dirProc '\Images\Intermediate_Products\glint_threshold_Per_Track.mat'],'Glint');

elseif isempty(whos('meanIm')) || isempty(whos('stdIm')) ||...
       isempty(whos('RM_Nr')) || isempty(whos('meanOriginal')) ||...
       isempty(whos('Glint'))
    
    % Load mean brightness and standard deviation imagesand glint threshold
    load([dirProc '\Images\Intermediate_Products\mean_Images_Per_Track.mat'])
    load([dirProc '\Images\Intermediate_Products\glint_threshold_Per_Track.mat'])

end

% Plot mean pixel brightness and its standard deviation for each track and 
% show examples of cropping high glint regions 
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

display_text('Step 5: Building and executing files for georeferencing images.','section')

if option_run_bat == 1

    %--- Build batch file for georeferencing images ---%

    % Loop through tracks
    for i=5 %1:length(tracks)
    
        % Check if track is stable
        if trackTag(i).stable==1
    
            % Compute mean northing and easting UTM coordinate of ith track and
            % convert it to decimal degree longitude and latitude
            xx=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1),'omitnan');
            yy=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2),'omitnan');
            [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
            
            % Set the paths to nongeoreferenced images and the output directory
            % for post-processed georeferenced images
            dirin_img = [dirProc 'Images\Intermediate_Products\nongeoreferenced_Images_With_Vignetting_Removed\Track_' num2str(i) '\']; 
            dirin_prj = dirProc;
            
            % Build batch file for running trimble georeferencing project 
            buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,i,res,option_plot,0);
        end
    end
    
    %--- Build batch file for georeferencing images for identifying high glint regions ---%
    
    % Loop through tracks 
    for i=5 %1:length(tracks)
    
        % Check if track is stable
        if trackTag(i).stable==1
            
            % Compute mean northing and easting UTM coordinate of ith track and
            % convert it to decimal degree longitude and latitude
            xx=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),1),'omitnan');
            yy=mean(A.data(tracks(i).Indices(1):tracks(i).Indices(2),2),'omitnan');
            [proj_lat,proj_lon]=utm2deg(xx,yy,utmZone);
            
            % Set the paths to masks for determing high glint and output
            % directory for georeferenced masks
            dirin_img = [dirProc 'Images\Intermediate_Products\mask_For_Determining_High_Glint_Areas\Track_' num2str(i) '\'];
            dirin_prj = dirProc;
            
            % Build batch file for running trimble georeferencing project
            buildBat(dirin_img,dirin_prj,proj_lat,proj_lon,i,res,option_plot,1);
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
end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine brightness threshold of breakers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_text('Step 6: Determining brightness threshold of breakers','section')

% Set directory to save georeferenced images
dirProc = 'D:\DEPLOYMENTS\PROGRAMS\WaveBreakingStats\data'; 
dirProcessed=[dirProc 'Output\'];

% Compute brightness threshold
[Threshold,N]=global_local_Thresh(peakPercentage,dirProcessed,tracks_Im,...
                                  trackTag,tracks,meanIm,RM_Nr,Glint,...
                                  localStep);

% Denoise local threshold estimate
if option_globalOrlocal == 0
    Threshold = denoiseThreshold(Threshold,tracks,trackTag);
end

% Save brightness threshold
save([dirProc 'Images\Intermediate_Products\brightness_threshold.mat'],'Threshold','N');

% Plot 
if option_plot == 1

    % Create a subdirectory for brightness threshold plots
    dirBT=[dirV '\Brightness_Threshold\'];
    if ~exist(dirBT, 'dir'), mkdir(dirBT); end

    % Generate plots
    plot_Thresh(Threshold,N,dirProcessed,dirV,D_Im,tracks_Im,trackTag,tracks,arrayHist,Glint,Glint_mask,meanIm,GlobalOrLocal);

end

display_text('Done!','body')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Lambda of c and compute wave breaking statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_text('Step 7: Determining Lambda of c and computing wave breaking statistics','section')

offSetEnds=0;

% The images are indexed in seconds from the begining of the week. This
% requires entering the start date manually.
StartTime=datenum([2021,5,2,0,0,0]);

% Determine wind direction from the data.
load('PLD2_1hz_ALL_1.mat');
WindTime=PLD2_1hz.time;
Speeds=PLD2_1hz.TWD;

% Loop through tracks
for i=1:length(tracks)
    
    % Check if flight track is stable 
    if trackTag(i).stable==1

        % Set 
        if i>=10
            dirV=[dirProc 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_' num2str(i) '\'];
        else
            dirV=[dirProc 'Lambdas_inbuilt_Threshold_0_pyramidScaled_V5\Track_0' num2str(i) '\'];
        end
        if ~exist(dirV, 'dir'), mkdir(dirV); end
        
       
        % 
        dirinIM=[dirProc 'Output\Track_' num2str(i) '\'];
        dirinIM_mask=[dirProc 'Output\TrackMask_' num2str(i) '\'];
        filenamesIM=dir([dirinIM '*.tif']);

        % 
        i1=1+3+(trackTag(i).range(1)+offSetEnds)-tracks(i).Indices(1);
        i2=1+(trackTag(i).range(2)-offSetEnds)-tracks(i).Indices(1)-1;

        % 
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(2);

        % 
        strN=filenamesIM(1).name;
        dxIm = str2num(strN(pozicija+1:pozicija+2))/100;
        
        % 
        WindDir=determine_wind(Speeds,WindTime,StartTime,tracks(i).Indices,An);

        % 
        process_Lambda_New(i1,i2-3,filenamesIM,dirinIM,dirV,Threshold(i).th,300,dxIm,0.2,WindDir,SpreadMax,Glint(i).list,dirinIM_mask,meanIm(i).im,RM_Nr(i).nr,GlobalOrLocal,trackTag(i),tracks(i))
    end
end


% Loop through tracks
for i=1:length(tracks)
    
    % Check if flight track is stable
    if trackTag(i).stable==1

        % 
        if i>=10
            dirL=[dirProc 'Lambdas_inbuilt\Track_' num2str(i) '\'];
        else
            dirL=[dirProc 'Lambdas_inbuilt\Track_0' num2str(i) '\'];
        end
        filenames=dir([dirL '*.mat']);
        
        % 
        dirinIM=[dirProc 'Output\Track_' num2str(i) '\'];
        filenamesIM=dir([dirinIM '*.tif']);
        
        % 
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(2);
        
        % 
        strN=filenamesIM(1).name;
        dxIm = str2num(strN(pozicija+1:pozicija+2))/100;

        % 
        pozicija=find('_'==filenamesIM(1).name,2);
        pozicija=pozicija(3);
        strN=filenamesIM(1).name;

        % 
        PARAM.dt = 0.2;
        PARAM.dx = str2num(strN(pozicija+1:pozicija+2))/100;
        
        % 
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

% Plot out other verification plots (wc coverage in images, total length of
% breaking, spectrogram of Lambda distribution).

%plotDists();

display_text('Done!','body')

%% Development Code

% Load buoy positions
% load('D:\MASS\Processed\L_Computations\L_Computations\Test\TFOEx21_DEP01_BuoyGPS.mat')

% Identify that EO files are processed through Trimble
% if isfolder([dirProc 'EO'])
%     movefile([dirProc 'EO'],[dirProc 'Trimble_EO']);
% end

% if option_globalOrlocal == 1
%     [Threshold,N]=global_Thresh(peakPercentage,dirProcessed,tracks_Im,trackTag,tracks,meanIm,RM_Nr,Glint);
% else 
%     [Threshold,N]=local_Thresh(peakPercentage,dirProcessed,tracks_Im,trackTag,tracks,meanIm,RM_Nr,Glint,localStep);
% end

%Threshold(i)=11000;
%process_Lambda_New(i1,i2,filenamesIM,dirinIM,dirV,Threshold(i),300,0.1,0.2)

%         dxIm = 0.3;

% % Plot out randomly chosen examples from the track
% dirV=[dirProc 'Quality_Control\LambdaIms_inbuilt\'];
% if ~exist(dirV, 'dir'), mkdir(dirV); end
% 
% dirV2=[dirProc 'Plots\DistributionsAll\'];
% if ~exist(dirV2, 'dir'), mkdir(dirV2); end
% 
% dirV2=[dirProc 'Plots\DistributionsAll_WC_Glint\'];
% if ~exist(dirV2, 'dir'), mkdir(dirV2); end
% 
% dirV2=[dirProc 'Plots\DistributionsAll_new\'];
% if ~exist(dirV2, 'dir'), mkdir(dirV2); end
% 
% dirV2=[dirProc 'Plots\Distributions_Inbuilt\'];
% if ~exist(dirV2, 'dir'), mkdir(dirV2); end
% 
% 
% dirV3=[dirProc 'Plots\Summary\'];
% if ~exist(dirV3, 'dir'), mkdir(dirV3); end
% 
% dirV4=[dirProc 'Plots\Alongtrack_OppositeDirection\'];
% if ~exist(dirV4, 'dir'), mkdir(dirV4); end
% 
% dirV4=[dirProc 'Plots\Alongtrack_test3\'];
% if ~exist(dirV4, 'dir'), mkdir(dirV4); end
% addpath('utilities')
% 
% 
% dirV2=[dirProc 'Plots\DistributionsAll_new\'];
% if ~exist(dirV2, 'dir'), mkdir(dirV2); end
% 

% % 
% windData='PLD2_1hz.mat';
% waveData='WG1mc_Snn_Hs_30min.mat';
% 
% % determine params for tfo
% for i=2:3%length(trackTag)
%     if trackTag(i).stable==1
%         
%         
%         WGTime=PLD2_AVG.time;
%         WGTime2=PLD2_FFT.time;
%         Hs=PLD2_AVG.Hs_fft;
%         Freq=0:0.01:0.82;
%         Spectra=[];
%         Spectra2=[];
%         for ii=1:83
%             Spectra(ii,:)=eval(['PLD2_FFT.avgSzz_' num2str(ii)]);
%             Spectra2(ii,:)=eval(['PLD2_FFT.avgSzz_' num2str(ii)]);
%         end
%         [Hs_all(i),freq_all(i)]=determine_params_TFO(Spectra,Freq,Hs,WGTime,WGTime2,StartTime,tracks(i).Indices,An,Spectra2);
%         
%         WindTime=PLD2_1hz.time;
%         us=movmean(PLD2_1hz.TWS*0.41/log(1.11/0.0002),100);
%         ustar(i)=determine_ustar(us,WindTime,StartTime,tracks(i).Indices,An);
%     end
% end
% 
% plotAll_TFO(dirV2,trackTag,StartTime,tracks,An,Hs_all,freq_all,ustar);
