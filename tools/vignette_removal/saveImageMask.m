function saveImageMask(dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirOut)
    
    %%%%
    % saveImageMask(dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,dirOut)
    %
    % Function for  
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
    %   meanIm    : Mean brightness of each pixel over the time period
    %               of the stable flight track (temporal mean). Mean
    %               brightness is smoothed using a 2D moving average.
    %   stdIm     : Sample standard deviation of brightness for each
    %               pixel over the time period of the stable flight 
    %               (temporal std). Standard deviation of brightness
    %               is smoothed using a 2D moving average.
    %   dirOut    : Path to directory for saving intermediate data products.
    %   sigma_ff  : Standard deviation of the Gaussian smoothing filter
    %               for the 2D image flate field correction and the 2D 
    %               Gaussian filtering of the squared deviations from the
    %               mean of the flat field image.
    %   n_sigma   : Number of standard deviations above or below the 
    %               median image brightness. Parameter is used for 
    %               determining which pixels are considered for the 
    %               calculation of the mean image brightness. 
    % 
    %   Returns
    %   -------
    %   Saves the equalized images in the directory specified by dirOut. 
    %  
    % 
    %%%%

    % This function crops out the areas with strong glint, removes vigneting,
    % and saves the processed images to dirProcessed folder.
    % Input: 
    % 

    numCores=8;
    poolobj = parpool(numCores);
    for i=14:17%10:4:14%1:length(tracks)
        if (trackTag(i).stable==1) & (~isempty(tracks(i).Indices))
            meanIm(i).im(isnan(meanIm(i).im))=0;
            stdIm(i).im(isnan(stdIm(i).im))=0;
            
            beginDif=tracks_Im(i).Indices(1);%+trackTag(i).range(1)-tracks(i).Indices(1);
            endDif=tracks_Im(i).Indices(2);%+trackTag(i).range(2)-tracks(i).Indices(2);
            medianIm=median(meanIm(i).im(:));
            medianStd=median(stdIm(i).im(:));
            
            dirV=[dirOut 'Track_' num2str(i) '\'];
    %         dirV=[dirProcessed 'Test\'];
            if ~exist(dirV, 'dir'), mkdir(dirV); end
            parfor j=beginDif:endDif
    %             a=imread([dirRaw D_Im(j).name]);
    %             a=double(im2gray(a));
    %             
    %             a2=a-meanIm(i).im+medianIm;
    %             a2=nanmean(a2(:))+(a2-nanmean(a2(:)))./stdIm(i).im*medianStd;
                %mina2(j+1-beginDif)=min(a2(:));
                % Shift by 5000 to eliminate any potential zeros.
                a2=uint16(meanIm(i).im/nanmax(meanIm(i).im(:))*256*256);
                
                imwrite(a2,[dirV D_Im(j).name]);
                
            end
        end
    end
    mina2=1;
    delete(poolobj)

end
