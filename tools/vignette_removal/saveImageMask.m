function saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirOut)
    
    %%%%
    % saveImageMask(D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirOut)
    %
    % Function for  
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
    %   stdIm     : Sample standard deviation of brightness for each
    %               pixel over the time period of the stable flight 
    %               (temporal std). Standard deviation of brightness
    %               is smoothed using a 2D moving average.
    %   dirOut    : Path to directory for saving intermediate data products.
    % 
    %   Returns
    %   -------
    %   Saves the image masks in the directory specified by dirOut. 
    %  
    % 
    %%%%

    % Initialize parallel computing if parallel computing is not already
    % running
    if isempty(gcp('nocreate'))

        % Set the number of CPU cores (i.e., the brains of the CPU that recieve
        % and execute operations for your computer) 
        numCores=8;
    
        % Run code on parallel pools for increases efficiency
        poolobj = parpool(numCores);

    else 

        % Initialize an empty poolobj
        poolobj = [];
    
    end

    % Loop through tracks
    for i=2 %1:length(tracks)

        % Check if track is stable and if full track has start and end time
        % indices
        if (trackTag(i).stable==1) && (~isempty(tracks(i).Indices))

            % Set nans in mean and standard deviation images to zero
            meanIm(i).im(isnan(meanIm(i).im))=0;
            stdIm(i).im(isnan(stdIm(i).im))=0;
            
            % Set beginning and end time indices for full flight track 
            beginDif=tracks_Im(i).Indices(1);                               %+trackTag(i).range(1)-tracks(i).Indices(1);
            endDif=tracks_Im(i).Indices(2);                                 %+trackTag(i).range(2)-tracks(i).Indices(2);
            
            % Create a subdirectory for image mask output after processing
            dirV=[dirOut 'Track_' num2str(i) '\'];
            if ~exist(dirV, 'dir'), mkdir(dirV); end

            % Loop through images 
            parfor j=beginDif:endDif
                
                % Normalize mean image by its maximum value and convert to
                % array elements to 16 bit unsigned integers  
                a2=uint16(meanIm(i).im/max(meanIm(i).im(:))*256*256);
                
                % Write processed image mask
                imwrite(a2,[dirV D_Im(j).name(1:end-4) '_Mask.tif']);
                
            end
        end
    end
    
    % End parallel computing session
    delete(poolobj)

end

%% Developmental Code
% This function crops out the areas with strong glint, removes vigneting,
% and saves the processed images to dirProcessed folder.
% Input: 
% 

%             % Compute the median mean and standard deviation values
%             medianIm=median(meanIm(i).im(:));
%             medianStd=median(stdIm(i).im(:));

%             a=imread([dirRaw D_Im(j).name]);
%             a=double(im2gray(a));
%             
%             a2=a-meanIm(i).im+medianIm;
%             a2=nanmean(a2(:))+(a2-nanmean(a2(:)))./stdIm(i).im*medianStd;
%mina2(j+1-beginDif)=min(a2(:));
% Shift by 5000 to eliminate any potential zeros.