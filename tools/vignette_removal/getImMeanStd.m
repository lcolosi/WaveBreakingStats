function [meanIm,stdIm,RM_Nr,meanOriginal] = getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,tracks,winSize,dirout)

    %%%%
    % [meanIm,stdIm,RM_Nr,meanOriginal]=getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,tracks,winSize,dirout)
    %
    % Function for  
    % 
    % First we obtain the mean of the image and of std for each track. This is
    % computed for the stable portion of the tracks, as indicated in trackTag
    % variable.
    % 
    %
    %   Parameters
    %   ----------
    %   dirRaw    : Path to raw (non-georeferenced) video images for a
    %               given flight.
    %   D_Im      : Filenames of the raw (non-georeferenced) video images.
    %   tracks_Im : Structure containing the start and end time indices of
    %               each full flight track (located in the indices field). 
    %               Indices derived from raw imagery file names. 
    %   trackTag  : Structure with three fields: 
    %               (1) stable: An identifier for whether the Track is 
    %                           stable or unstable based on the inputted
    %                           roll, pitch, and heading criteria. Here, we
    %                           denote stability with the following
    %                           convention:
    %                                   0 -> Unstable
    %                                   1 -> Stable    
    %               (2) range: The start and end time indices of the stable
    %                          flight track. 
    %               (3) stats: An array containing the statistics from the
    %                          initial stability analysis using the maxPer 
    %                          and shift parameters. The array has the
    %                          following structure: 
    %                                  [mean_r, std_r;... 
    %                                   mean_p, std_p;...
    %                                   mean_h, std_h; 
    %                                   rmN(1), rmN(2)]
    %                           where mean and standard deviations are
    %                           computed over the interval: 
    %                                   tracks(i).Indices(1)+rmN(1) to 
    %                                   tracks(i).Indices(2)-rmN(2) 
    %   tracks    : Structure containing the start and end time indices of
    %               each full flight track (located in the indices field).
    %               Indices derived from EO file.
    %   winSize   : 
    %   dirout    : 
    % 
    %   Returns
    %   -------
    %   meanIm       : 
    %   stdIm        : 
    %   RM_Nr        :
    %   meanOriginal : 
    % 
    % 
    % The RM_Nr and meanOriginal variables check for significant variations in
    % image lighting. RM_Nr represents image numbers to be removed, and
    % meanOriginal is the mean lighting in every image
    % 
    % This functions requires the Parallel computing toolbox to access the
    % parpool function.  
    % 
    % Version: 5.0 from Teodor's Codes for Lambda in TFO_2021 directory in 
    % airseaserver28 
    % 
    %%%%

    % Set the number of CPU cores (i.e., the brains of the CPU that recieve
    % and execute operations for your computer) 
    %numCores=6;

    % Run code on parallel pools for increases efficiency
    %poolobj = parpool(numCores);

    % Add utilities directory to the top of the current working directory 
    %addpath('utilities');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute mean brightness for every pixel on the track
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: This is only done if a tracks is identified as stable

    % Loop through tracks 
    for i=1:length(tracks_Im)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Remove vignetting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove vignetting, cut sunglint using the previous iteration,
        % determine brightness threshold for breaking.

        % Check if track is stable 
        if trackTag(i).stable==1

            % Set beginning and end indices for stable flight period 
            beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
            endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);
            
            %--- Check image format (8 or 16 bit) using first image file---%

            % Load first raw image from ith track's stable period
            a=imread([dirRaw D_Im(beginDif).name]);
            
            % Set pixel value range for 8 and 16 bit images
            if isa(a,'uint16')
                arrayHist=0:16:256*256;
            elseif isa(a,'uint8')
                arrayHist=0:256;
            else
                % Print warning for unknown image format
                fprintf('Unknown image format. Only 8 or 16 bit images are recognized');
            end

            % 
            [meanIM,stdIM,meanOriginal(i).nr]=detsPar(beginDif,endDif,dirRaw,D_Im);
            [tempp,RM_Nr(i).nr]=rmoutliers(meanOriginal(i).nr(beginDif:endDif));
            [N,Threshold]=detsPar2(beginDif,endDif,dirRaw,D_Im,meanIM,stdIM,arrayHist,RM_Nr(i).nr);
            
            
            [meanIm(i).im,stdIm(i).im]=detsPar3(beginDif,endDif,dirRaw,D_Im,meanIM,stdIM,arrayHist,Threshold,RM_Nr(i).nr);
            
            [meanOriginal(i).nr,RM_Nr(i).nr,meanIm(i).im,stdIm(i).im]=determineMeanStd(beginDif,endDif,dirRaw,D_Im,arrayHist);
            
            meanIm(i).im=mean2D(meanIm(i).im,winSize);
            stdIm(i).im=mean2D(stdIm(i).im,winSize);
        
        % If the start and end time indices for the ith track are empty,
        % set mean and standard deviation of image to empty arrays
        else
            meanIm(i).im=[];
            stdIm(i).im=[];
        end
    end
    
    % End parallel computing session
    delete(poolobj)
    return

%% Development Code
% This variable tracks if image format is recognised (8 or 16 bit)
% Err=0;
%     if Err==1
%         break
%     end
% 

%     else
%         fprintf('Unknown image format. 8/16 bit are recognized');
%         Err=1;
%         %break
%     end

% Convert RGB image to grayscale and pixel values to double
% a=double(im2gray(a));