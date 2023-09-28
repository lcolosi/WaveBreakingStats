function [meanOriginal,RM_Nr,meanIm,stdIm] = determineMeanStd(beginDif,endDif,dirRaw,D_Im,sigma_ff,B_threshold,n_sigma)

    %%%%
    % [meanOriginal,RM_Nr,meanIm,stdIm] = determineMeanStd(beginDif,endDif,dirRaw,D_Im,arrayHist,sigma_ff,B_threshold,n_sigma)
    %
    % Function for determining the mean and standard deviation of 
    % brightness for each pixel over the time period of the stable flight
    % track. Outlier images (i.e., images with significantly high or low 
    % mean brightness) are excluded from calculations.  
    %
    %   Parameters
    %   ----------
    %   beginDif    : Beginning time indices for stable flight period. 
    %   endDif      : End time indices for stable flight period. 
    %   dirRaw      : Path to raw (non-georeferenced) video images for a
    %                 given flight.
    %   D_Im        : Filenames of the raw (non-georeferenced) video images.
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
    %   meanIm       : Mean brightness of each pixel over the time period
    %                  of the stable flight track (temporal mean). Here, at
    %                  each time step (that are not identified as a 
    %                  brightness outlier by by RM_Nr), pixels greater than
    %                  n_sigma are not included in the calculation of the 
    %                  mean. 
    %   stdIm        : Sample standard deviation of brightness for each
    %                  pixel over the time period of the stable flight 
    %                  (temporal std). Same calculation criteria for the
    %                  mean brightness is applied to the standard
    %                  deviation calculation.  
    %   RM_Nr        : A logical array indicating the images whose lighting is
    %                  significantly higher or lower than the median 
    %                  brightness of all images (i.e. the outliers). 
    %                  Significance deviation from the median is defined as 
    %                  more than three scaled median absolute deviations 
    %                  from the median  This is the quantity that identifies 
    %                  significant variations in an image's lighting.
    %   meanOriginal : Mean brightness over all pixels for each image 
    %                  (spatial average). Here, the mean of the all pixels
    %                  below the brightness threshold set by B_threshold of
    %                  each image is computed. This is computed to check
    %                  for significant variations in image lighting.
    % 
    %   Notes
    %   -----
    %   This functions requires the Parallel computing and Image processing 
    %   toolbox for access to the parfor and imflatfield functions
    %   respectively.
    % 
    %   Important: You can not use the debug mode in parfor loops! 
    % 
    %%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute mean brightness over all pixels for each image 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate waitbar
    pw = PoolWaitbar(endDif-beginDif, 'Computing mean brightness over all pixels for each image.');

    % Loop through images 
    parfor j=beginDif:endDif
        
        % Update waitbar
        increment(pw)

        % Load jth image
        a=imread([dirRaw D_Im(j).name]);
        
        % Convert RGB image to grayscale and pixel values to double
        a=double(im2gray(a));
    
        % Sort image brightness elements in ascending order
        ImSeries2=sort(a(:));
        
        % Remove the highest 20% of pixels from mean brightness
        % computation. 
        ImSeries2=ImSeries2(1:floor(B_threshold*length(ImSeries2)));
        
        % Compute the mean brightness of the jth image (mean of pixels)
        meanOriginal(j)=mean(ImSeries2);

    end

    % Obtain a logical array indicating the images whose lighting is
    % significantly higher or lower than the median brightness of all
    % images (i.e. the outliers). 
    [~,RM_Nr]=rmoutliers(meanOriginal(beginDif:endDif));

    % Obtain indices of outliers 
    idx = beginDif:endDif;               
    idx_outliers = idx(RM_Nr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the mean brightness of each pixel over the stable track period
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Read first image from stable flight period
    aTemp=imread([dirRaw D_Im(beginDif).name]);

    % Initialize zeros arrays for mean image and counter
    meanIM2=zeros(size(aTemp));
    counter=zeros(size(aTemp));

    % Generate another waitbar
    pw = PoolWaitbar(endDif-beginDif, 'Computing the mean brightness of each pixel over the stable track period.');
    
    % Loop through images
    parfor j=beginDif:endDif

        % Check if jth image is identified as a brightness outlier
        if ~ismember(j,idx_outliers)

            % Update waitbar
            increment(pw)
            
            % Load jth image
            a=imread([dirRaw D_Im(j).name]);

            % Convert RGB image to grayscale and pixel values to double
            a=double(im2gray(a));

            % Preform a flat-field correction (a digital image technique 
            % to mitigate the image detector pixel-to-pixel sensitivity and
            % distortions in the optical path)
            a2=imflatfield(a,sigma_ff);
            
            % Set pixels more than n_sigma standard deviations above or below the
            % median brightness of the image to NaN (not considered in the 
            % mean image brightness for the track)
            a2(a2>median(a2(:))+n_sigma*std(a2(:)) | a2<median(a2(:))-n_sigma*std(a2(:)))=nan;

            % Add to counter at pixels within n_sigma standard deviations
            % from the median
            counter=counter+double(~isnan(a2));
            
            % Set all pixels outside n_sigma standard deviations from the
            % mean to zero in the raw image
            a(isnan(a2))=0;
            
            % Sum the jth image
            meanIM2=meanIM2+a;
        end
    end

    % Compute the mean image brightness for the ith track
    meanIm=meanIM2./counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the standard deviation of brightness for each pixel over the stable track period
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize zeros arrays for mean image and counter 
    stdIM2=zeros(size(aTemp));
    counter=zeros(size(aTemp));
    
    % Generate another waitbar
    pw = PoolWaitbar(endDif-beginDif, 'Computing the standard deviation of brightness for each pixel over the stable track period.');
    
    % Loop through images
    parfor j=beginDif:endDif

        % Check if jth image is identified as a brightness outlier
        if ~ismember(j,idx_outliers)                                        % Old code: if ~ismember(j,RM_Nr-1+beginDif)

            % Update waitbar
            increment(pw)

            % Load jth image
            a=imread([dirRaw D_Im(j).name]);
            
            % Convert RGB image to grayscale and pixel values to double
            a=double(im2gray(a));

            % Preform a flat-field correction
            a2=imflatfield(a,sigma_ff);
    
            % Set pixels more than n_sigma standard deviations above or below the
            % median brightness of the image to NaN (not considered in the 
            % mean image brightness for the track)
            a2(a2>median(a2(:))+n_sigma*std(a2(:)) | a2<median(a2(:))-n_sigma*std(a2(:)))=nan;
            
            % Add to counter at pixels within n_sigma standard deviations
            % from the median
            counter=counter+double(~isnan(a2));
            
            % Compute the square of the image minus the mean track image
            addTemp=((a-meanIm).^2);
            
            % Set all pixels outside n_sigma standard deviations from the
            % mean to zero in the squared difference quantity 
            addTemp(isnan(a2))=0;
            
            % Sum the squared differences
            stdIM2=stdIM2+addTemp;
        end
    end

    % Compute sample standard deviation
    stdIm=(stdIM2./(counter-1)).^0.5;                                       % Old code: stdIm=(stdIM2./counter).^0.5;
end
