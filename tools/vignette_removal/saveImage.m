function saveImage(dirRaw,D_Im,tracks_Im,trackTag,tracks,dirout,sigma_ff,n_sigma,meanIm,stdIm,meanOriginal,sigBrightness)

    %%%%
    % saveImage(dirRaw,D_Im,tracks_Im,trackTag,tracks,dirout,sigma_ff,n_sigma,meanIm,stdIm,meanOriginal,sigBrightness)
    %
    % Function for removing vingetting and equalizing images during the 
    % stable flight period. Vignetting is removed using a flat field 
    % correction. The image is further equalized by normalizing each
    % pixel's deviation from the mean brightness by the 2D gaussian filter
    % of the squared deviations from the mean (this normalization quantity
    % is itself normalized by its mean value). The final equalized image is
    % written to a new graphics file with the same name as the raw image. 
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
    %   dirout    : Path to directory for saving intermediate data products.
    %   sigma_ff  : Standard deviation of the Gaussian smoothing filter
    %               for the 2D image flate field correction and the 2D 
    %               Gaussian filtering of the squared deviations from the
    %               mean of the flat field image. (originally 450)
    %   n_sigma   : Number of standard deviations above or below the 
    %               median image brightness. Parameter is used for 
    %               determining which pixels are considered for the 
    %               calculation of the mean image brightness.
    %   meanIm    : Mean brightness of each pixel over the time period
    %               of the stable flight track (temporal mean). Mean
    %               brightness is smoothed using a 2D moving average.
    %   stdIm     : Sample standard deviation of brightness for each
    %               pixel over the time period of the stable flight 
    %               (temporal std). Standard deviation of brightness
    %               is smoothed using a 2D moving average.
    %   meanOriginal  : Mean brightness over all pixels for each image 
    %                   (spatial average). Here, the mean of the all pixels
    %                   below the brightness threshold set by B_threshold of
    %                   each image is computed. This is computed to check
    %                   for significant variations in image lighting. 
    %   sigBrightness : Maximum allowed standard deviation of meanOriginal 
    %                   to constitute constant brightness along the track. 
    %                   For tracks with std(meanOriginal) >= sigBrightness, 
    %                   the variable lighting algorithm for removing
    %                   vignetting/equalizing image is used. Otherwise, the
    %                   constant lighting algorithm is used. 
    % 
    %   Returns
    %   -------
    %   Saves the equalized images in the directory specified by dirout. 
    %  
    % 
    %%%%

    % Initialize parallel computing if parallel computing is not already
    % running
    if isempty(gcp('nocreate'))

        % Set the number of CPU cores (i.e., the brains of the CPU that recieve
        % and execute operations for your computer) 
        numCores=10; %feature('numcores');
    
        % Run code on parallel pools for increases efficiency
        poolobj = parpool(numCores);

    else

        % Initialize an empty poolobj
        poolobj = gcp('nocreate');

    end

    % Generate waitbar
    pos = [585 421.8750 270 56.2500];
    f = waitbar(0,'Please wait...','Position', [pos(1) pos(2)+2*pos(4) pos(3) pos(4)]);

    % Loop through tracks 
    for i=1 %1:length(tracks)

        % Update waitbar
        waitbar(i/length(tracks),f,...
            {['Running saveImage.m : On track ' num2str(i) ' of ' num2str(length(tracks))]; [num2str(round((i/length(tracks))*100)) '$\%$ complete...']})

        % Check if track is stable and if full track has start and end time
        % indices
        if (trackTag(i).stable==1) && (~isempty(tracks(i).Indices))

            % Set beginning and end indices for stable flight period 
            beginSP = trackTag(i).range(1);
            endSP = trackTag(i).range(2);

            % Set the mean and standard deviation pixel images to avoid 
            % high data communication overhead in parfor loop 
            meanIm_temp = meanIm(i).im;
            stdIm_temp = stdIm(i).im;

            % Set nans in mean and standard deviation images to zero
            meanIm_temp(isnan(meanIm_temp))=0;
            stdIm_temp(isnan(stdIm_temp))=0;
            
            % Compute the median mean and standard deviation values
            medianIm=median(meanIm_temp(:));
            medianStd=median(stdIm_temp(:));
            
            % Create a subdirectory for image output after processing
            dirV=[dirout 'Track_' num2str(i) '\'];
            if ~exist(dirV, 'dir'), mkdir(dirV); end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Remove vignetting and equalize image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Generate waitbar
            pw = PoolWaitbar(endSP-beginSP, 'Removing vignetting and equalizing image.');

            % Loop through images
            parfor j=beginSP:endSP
                
                % Update waitbar
                increment(pw)

                % Load jth image
                a=imread([dirRaw D_Im(j).name]);
                
                % Convert RGB image to grayscale and pixel values to double
                a=double(im2gray(a));

                %--- Approach 1: Variable lighting conditions along track ---%
                if std(meanOriginal(i).nr) >= sigBrightness                 %#ok

                    % Preform a flat-field correction (helps remove vignetting)
                    a2=imflatfield(a,sigma_ff);
    
                    % Compute mean brightness over all pixels
                    meanVal = mean(a2(:),'omitnan');
                    
                    % Copy a2 variable to a new variable a3
                    a3=a2;
                    
                    % Set pixels more than n_sigma standard deviations above or
                    % below the median brightness of the image to image's mean
                    % brightness value 
                    a3(a3>median(a3(:))+n_sigma*std(a3(:)) | a3<median(a3(:))-n_sigma*std(a3(:)))=meanVal;
    
                    % Set size of the Gaussian filter 
                    filterSize = 2*(2*sigma_ff)+1;
    
                    % Compute the 2D Gaussian filtering of squared deviations
                    % from the mean. 
                    squareIm = imgaussfilt((a3-meanVal).^2, sigma_ff, 'Padding', 'symmetric', 'FilterSize', filterSize); 
                    
                    % Compute the population variance like quantity 
                    % for the gaussian filtered squared deviations for
                    % normalization of 2D Gaussian filter of square deviations
                    % from the mean 
                    meanVal2 = mean(squareIm(:),'omitnan');
    
                    % Normalize the deviations from the mean by the 2D gaussian
                    % filtered squared deviations from the mean of the image 
                    % in order to equalize image. 
                    a2 = meanVal+(a2-meanVal)./sqrt(squareIm/meanVal2);

                %--- Approach 2: Constant lighting conditions along track ---%
                else

                    % Equalize image by removing mean pixel brightness and
                    % replacing it with the median value of the mean
                    % brightness image
                    a2 = a-meanIm_temp+medianIm;
            
                    % Further equalize image by normalizing the deviation
                    % from the mean by the standard deviation of pixel
                    % brightness and replacing it with median standard
                    % deviation pixel brightness value
                    a2 = mean(a2(:),'omitnan')+((a2-mean(a2(:),'omitnan'))./stdIm_temp)*medianStd;

                end
                
                % Shift by 5000 to eliminate any potential zeros (tif files
                % cannot contain any zero or negative valued elements)
                a2=uint16(a2+5000);
                
                % Write processed image 
                imwrite(a2,[dirV D_Im(j).name]);
                
            end
        end
    end

    % End parallel computing session and close wait bar
    delete(poolobj); close(f)

end
