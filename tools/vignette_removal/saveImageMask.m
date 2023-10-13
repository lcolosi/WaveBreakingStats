function saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,dirOut)
    
    %%%%
    % saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,dirOut)
    %
    % Function for identifying high sun glint regions in each image during 
    % the stable flight period for each track. More specifically, this function 
    % normalizes each image for a given track by the maximum brightness
    % in the image and saves them in a new tif file. The output of this
    % function will be georeferenced using trimble. Then using the 
    % glint brightness threshold from det_glint.m script, the high sun
    % glint region will be identified for each georeferenced image. The
    % mask associated with this region will be applied to the
    % georeferenced image with vignetting removed.      
    %
    %   Parameters
    %   ---------- 
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
    %   meanIm    : Mean brightness of each pixel over the time period
    %               of the stable flight track (temporal mean). Mean
    %               brightness is smoothed using a 2D moving average.
    %   dirOut    : Path to directory for saving intermediate data products.
    % 
    %   Returns
    %   -------
    %   Saves the normalized images in the directory specified by dirOut. 
    %  
    % 
    %%%%

    % Initialize parallel computing if parallel computing is not already
    % running
    if isempty(gcp('nocreate'))

        % Set the number of CPU cores (i.e., the brains of the CPU that recieve
        % and execute operations for your computer) 
        numCores=feature('numcores');
    
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
    for i=5 %1:length(tracks)

        % Update waitbar
        waitbar(i/length(tracks),f,...
            {['Running saveImageMask.m : On track ' num2str(i) ' of ' num2str(length(tracks))]; [num2str(round((i/length(tracks))*100)) '$\%$ complete...']})

        % Check if track is stable and if full track has start and end time
        % indices
        if (trackTag(i).stable==1) && (~isempty(tracks(i).Indices))

            % Set the mean pixel images to avoid high data communication
            % overhead in parfor loop 
            meanIm_temp = meanIm(i).im;

            % Set nans in mean and standard deviation images to zero
            meanIm_temp(isnan(meanIm_temp))=0;
            
            % Set beginning and end time indices for the stable flight period 
            beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
            endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);
            
            % Create a subdirectory for image mask output after processing
            dirV=[dirOut 'Track_' num2str(i) '\'];
            if ~exist(dirV, 'dir'), mkdir(dirV); end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Generate mask for identifying high sun-glint regions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Generate waitbar
            pw = PoolWaitbar(endDif-beginDif, 'Generating mask for identifying high sun-glint regions.');

            % Loop through images 
            parfor j=beginDif:endDif
                
                % Update waitbar
                increment(pw)

                % Normalize mean image by its maximum value and convert to
                % array elements to 16 bit unsigned integers  
                a2=uint16((meanIm_temp/max(meanIm_temp(:)))*256*256);
                
                % Write processed image mask
                imwrite(a2,[dirV D_Im(j).name(1:end-4) '_Mask.tif']);
                
            end
        end
    end
    
    % End parallel computing session
    delete(poolobj)

end
