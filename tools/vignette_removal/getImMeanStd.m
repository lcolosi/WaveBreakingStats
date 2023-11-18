function [meanIm,stdIm,RM_Nr,outliers,meanOriginal] = getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,winSize,sigma_ff,B_threshold,n_sigma)

    %%%%
    % [meanIm,stdIm,RM_Nr,outliers,meanOriginal]=getImMeanStd(dirRaw,D_Im,tracks_Im,trackTag,winSize,sigma_ff,B_threshold,n_sigma)
    %
    % Function for computing the mean and standard deviation of each pixel
    % for each track during the stable periods. Additionally outliers are 
    % identified in the mean image. 
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
    %   winSize   : Window or Kernel size for 2D move average. Window is
    %               square so winSize is a scale quantity representing the
    %               side length of the square window.  
    %   sigma_ff    : Standard deviation of the Gaussian smoothing filter
    %                 for the 2D image flate field correction. 
    %   B_threshold : Fraction of the pixels (sorted in ascending 
    %                 numerical order according to brightness) that are 
    %                 considered in the mean pixel brightness calculation. 
    %                 For example, if per_thresh = 0.8, then the top 20% of
    %                 pixel values are excluded from the mean pixel
    %                 calculation.
    %   n_sigma     : Number of standard deviations above or below the 
    %                 median image brightness. Parameter is used for 
    %                 determining which pixels are considered for the 
    %                 calculation of the mean image brightness.
    % 
    %   Returns
    %   -------
    %   meanIm       : Mean brightness of each pixel over
    %                  the time period of the stable flight track 
    %                  (temporal mean). Mean brightness is smoothed using
    %                  a 2D moving average. 
    %   stdIm        : Sample standard deviation of brightness for each
    %                  pixel over the time period of the stable flight 
    %                  (temporal std). Standard deviation of brightness
    %                  is smoothed using a 2D moving average.  
    %   RM_Nr        : A logical array indicating the images whose lighting is
    %                  significantly higher or lower than the median 
    %                  brightness of all images (i.e. the outliers). 
    %                  Significance deviation from the median is defined as 
    %                  more than three scaled median absolute deviations 
    %                  from the median  This is the quantity that identifies 
    %                  significant variations in an image's lighting.
    %   outliers     : A structure containing three fields:
    %                       (1) L : The lower threshold value for outliers. 
    %                               This value is defined as
    %                               three scaled median absolute 
    %                               deviations below the median.  
    %                       (2) U : The upper threshold value for outliers. 
    %                               This value is defined as
    %                               three scaled median absolute 
    %                               deviations above the median.
    %                       (3) C : The median of mean image brightness
    %                               time series defined here as the center
    %                               value of the outlier analysis.
    %                       (4) N : The total number of images identified
    %                               as having significant deviation in
    %                               mean image brightness.
    %   meanOriginal : Mean brightness over all pixels for each image 
    %                  (spatial average). Here, the mean of the all pixels
    %                  below the brightness threshold set by B_threshold of
    %                  each image is computed. This is computed to check
    %                  for significant variations in image lighting. 
    %  
    %   Notes
    %   -----
    %   This functions requires the Parallel computing toolbox to access the
    %   parpool function. This is only accessible through the air-sea lab's
    %   matlab account.   
    % 
    %   Version: 5.0 from Teodor's Codes for Lambda in TFO_2021 directory in 
    %   airseaserver28 
    % 
    %%%%
    
    % Initialize parallel computing if parallel computing is not already
    % running
    if isempty(gcp('nocreate'))                                             % || gcp('nocreate').NumWorkers < feature('numcores')

        % Set the number of CPU cores (i.e., the brains of the CPU that recieve
        % and execute operations for your computer) 
        numCores=10; %feature('numcores');
    
        % Run code on parallel pools for increases efficiency
        poolobj = parpool(numCores);

    else 

        % Initialize an empty poolobj
        poolobj = gcp('nocreate');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute mean and std of brightness for each pixel over stable track period
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: This is only done if a tracks is identified as stable

    % Generate waitbar
    pos = [585 421.8750 270 56.2500];
    f = waitbar(0,'Please wait...','Position', [pos(1) pos(2)+2*pos(4) pos(3) pos(4)]);

    % Loop through tracks 
    for i=1 %1:length(tracks_Im)

        % Update waitbar
        waitbar(i/length(tracks_Im),f,...
            {['Running getImMeanStd.m : On track ' num2str(i) ' of ' num2str(length(tracks_Im))]; [num2str(round((i/length(tracks_Im))*100)) '$\%$ complete...']})

        % Check if track is stable 
        if trackTag(i).stable==1

            % Set beginning and end indices for stable flight period 
            beginSP = trackTag(i).range(1);
            endSP = trackTag(i).range(2);

            % Compute the mean and standard deviation of the image plus identify outliers       
            [meanOriginal(i).nr,RM_Nr(i).nr,outliers(i).nr, meanIm(i).im,stdIm(i).im]=determineMeanStd(beginSP,endSP,dirRaw,D_Im,sigma_ff,B_threshold,n_sigma); %#ok
            
            % Compute the 2D moving average of the mean and standard deviation of brightness 
            meanIm(i).im=mean2D(meanIm(i).im,winSize);                      %#ok
            stdIm(i).im=mean2D(stdIm(i).im,winSize);                        %#ok
        
        % If the start and end time indices for the ith track are empty,
        % set mean and standard deviation of image to empty arrays
        else
            meanIm(i).im=[];                                                %#ok
            stdIm(i).im=[];                                                 %#ok
        end
    end

    % Close waitbar
    close(f)
    
    % End parallel computing session
    delete(poolobj)

end
