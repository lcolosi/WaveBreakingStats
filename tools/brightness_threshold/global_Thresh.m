function [Threshold,N]=global_Thresh(peakPercentage,dirProcessed,tracks_Im,trackTag,tracks,meanIm,RM_Nr,Glint)

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
            {['Running global_Thresh.m : On track ' num2str(i) ' of ' num2str(length(tracks))]; [num2str(round((i/length(tracks))*100)) '$\%$ complete...']})
        
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

            % Initialize empty array for 
            N(i).hi=[];

            % Generate waitbar
            pw = PoolWaitbar(endDif-beginDif, 'Removing vignetting and equalizing image.');
            
            % Loop through images
            parfor j=beginTrack:endTrack
                
                % Check if image has significant variation in lighting 
                % (do not process if significant variations are present)
                if RM_Nr(i).nr(j-beginTrack-1)==0
 
                    % Update waitbar
                    waitbar(j/(endTrack-beginTrack),f_image,...
                        {''; [num2str(round((i/endTrack-beginTrack)*100)) '$\%$ complete...']})

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


                    % 
                    if isempty(N(i).hi)
                        N(i).hi=histcounts(a,arrayHist);
                    else
                        N(i).hi=N(i).hi+histcounts(a,arrayHist);
                    end
                    
                    
                    T_temp(j,:)=histcounts(a,arrayHist);
                    
                end
            end
            

            Threshold(i).th=det_Threshold(N(i).hi,arrayHist,peakPercentage);
 
        end
    end
end

%% Developmental code

%     imTemp0=mean2D(meanIm(i).im,15);
%     imTemp=sort(imTemp0(:));
%     imTemp=imTemp(1:floor(end/2));
%     imMedian=median(imTemp);
%     imStd=std(imTemp);



% Find biggest connected area corresponding to glint using the
% selected criteria.
%     BinaryImage=imTemp0>imMedian+imStd*3;
%     stats = regionprops(BinaryImage, 'Area','PixelIdxList');
%     [biggest_area,Indeks] = max( [stats.Area] );
%     Glint(i).im=stats(Indeks).PixelIdxList;

%a=imread([dirProcessed '\Track_' num2str(i) '\' D_Im(j).name]);

%assignin('base','imageTemp',a_mask);
%assignin('base','imageTemp',a);
%assignin('base','test',Threshold);

% beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
% endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);