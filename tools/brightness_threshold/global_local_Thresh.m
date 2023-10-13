function [Threshold,N] = global_local_Thresh(peakPercentage,dirProcessed,tracks_Im,trackTag,tracks,meanIm,RM_Nr,Glint)

    %%%%
    % [Threshold,N]=global_Thresh(peakPercentage,dirProcessed,tracks_Im,trackTag,tracks,meanIm,RM_Nr,Glint)
    %
    % Function for determining the brightness theshold (over which a pixel is 
    % considered to be whitecap).    
    %
    %   Parameters
    %   ---------- 
    %   peakPercentage  : 
    %   dirProcessed
    % 
    %   Returns
    %   ------- 
    %    
    %%%%

    % Initialize parallel computing if parallel computing is not already
    % running
    if isempty(gcp('nocreate'))

        % Set the number of CPU cores (i.e., the brains of the CPU that recieve
        % and execute operations for your computer) 
        numCores= 10; %feature('numcores');
    
        % Run code on parallel pools for increases efficiency
        poolobj = parpool(numCores);

    else

        % Initialize an empty poolobj
        poolobj = gcp('nocreate');

    end

    % Set pixel value range for 16 bit images
    arrayHist=0:16:256*256;

    % Generate waitbar
    pos = [585 421.8750 270 56.2500];
    f = waitbar(0,'Please wait...','Position', [pos(1) pos(2)+2*pos(4) pos(3) pos(4)]);

    % Loop through tracks
    for i=5 %1:length(tracks)

        % Update waitbar
        waitbar(i/length(tracks),f,...
            {['Running global_local_Thresh.m : On track ' num2str(i) ' of ' num2str(length(tracks))]; [num2str(round((i/length(tracks))*100)) '$\%$ complete...']})
        
        % Check if track is stable
        if trackTag(i).stable==1
            
            % Set path to georeferenced images and obtain file names
            dirP=[dirProcessed 'Track_' num2str(i) '\'];
            filenames=dir([dirP '*.tif']);
            
            % Set path to georeferenced image for high sun-glint masking
            % and obtain file names 
            dirP_mask=[dirProcessed 'TrackMask_' num2str(i) '\'];
            filenames_mask=dir([dirP_mask '*.tif']);

            % Set beginning and end indices for stable flight period
            beginTrack = trackTag(i).range(1)-tracks(i).Indices(1) + 1;
            endTrack = length(filenames)+trackTag(i).range(2)-tracks(i).Indices(2); 

            % Initialize an empty array for image histograms
            N_image = [];

            % Generate waitbar
            pw = PoolWaitbar(endDif-beginDif, 'Removing vignetting and equalizing image.');
            
            % Loop through images
            parfor j=beginTrack:endTrack
                
                % Check if image has significant variation in lighting 
                % (do not process if significant variations are present)
                if RM_Nr(i).nr(j-beginTrack-1)==0
 
                    % Update waitbar
                    increment(pw)

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Remove regions of high sun-glint
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Load jth images
                    a = imread([dirP filenames(j).name]);
                    a_mask = imread([dirP_mask filenames_mask(j).name]);

                    % Convert RGB images to grayscale and pixel values to double
                    a=double(im2gray(a));
                    a_mask=double(im2gray(a_mask));

                    % 
                    a_mask=a_mask*max(meanIm(i).im(:),[],'omitnan')/256/256;

                    % Obtain glint magnitude threshold
                    glintMag=Glint(i).list;

                    % Generate a logical array for identifying pixels with brightness
                    % above the glint threshold
                    BinaryImage=a_mask > glintMag;

                    % Calculate the following properties for each region
                    % with brightness greater than the glint threshold: 
                    %   (1) Area: Number of pixels in the region (returns a
                    %             scalar for each region identified).
                    %   (2) PixelIdxList : Linear indices of the pixels in the
                    %                     region (returns a vector of p 
                    %                     elements long for each region where
                    %                     p is the total number of pixels in
                    %                     the region)
                    % The statss variable is a Nx1 structure with fields Area
                    % and PixelIdxList where N is the number of regions
                    % identified.
                    statss = regionprops(BinaryImage, 'Area','PixelIdxList');

                    % Find the region with the largest area
                    [~,Indeks] = max( [statss.Area] );

                    % Obtain the pixel indices of this region; this constitutes
                    % the high glint mask
                    GlintTemp=statss(Indeks).PixelIdxList;

                    % Apply glint mask to image
                    a(GlintTemp)=nan;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Compute histograms for each image
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    
                    % Counts number of pixel brightness values within bins
                    % whose edges are specified by arrayHist for jth image
                    N_image(j,:)=histcounts(a,arrayHist);
                    
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Determine brightness threshold of breakers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            %--- Global ---% 
            if option_globalOrlocal == 1

                % Sum over all histograms 
                N(i).hi = sum(N_image,1); 

                % Compute global brightness threshold 
                Threshold(i).th = det_Threshold(N(i).hi,arrayHist,peakPercentage);

            %--- Local ---% 
            else 

                % 
                N(i).hi = movmean(N_image,localStep,1,'omitnan');

                % Loop through images
                for j=beginTrack:endTrack
    
                    % Check if image has significant variation in lighting 
                    % (do not process if significant variations are present)
                    if RM_Nr(i).nr(j-initalNr)==0
    
                        % 
                        Threshold(i).th(j)=det_Threshold(N(i).hi(j,:),arrayHist,peakPercentage);
                    end
                end
            end
        end
    end
end

%% Development code 

% % 
% if isempty(N(i).hi)
%     N(i).hi=histcounts(a,arrayHist);
% else
%     N(i).hi=N(i).hi+histcounts(a,arrayHist);
% end
