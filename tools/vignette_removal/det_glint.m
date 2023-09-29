function [Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag)

    %%%%
    % [Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag)
    %
    % Function for determining the glint brightness threshold and the
    % corresponding high-glint mask for the track.   
    %
    %   Parameters
    %   ---------- 
    %   meanIm   : Mean brightness of each pixel over the time period
    %              of the stable flight track (temporal mean). Mean
    %              brightness is smoothed using a 2D moving average.
    %   stdMag   : 
    %   trackTag : Structure with three fields: 
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
    %   tracks   : Structure containing the start and end time indices of
    %              each full flight track (located in the indices field).
    %              Indices derived from EO file.
    % 
    %   Returns
    %   -------
    %   Glint      : Glint threshold brightness for each track. A scalar
    %                quantity defined as the median value of mean pixel
    %                brightness image plus n times the median value of the 
    %                standard deviation of pixel brightness.     
    % 
    %%%%

    % Loop through tracks
    for i=1 %1:length(tracks)
        
        % Check if track is stable 
        if trackTag(i).stable==1
            
            % Sort mean pixel brightness in ascending numerical order
            imTemp=sort(meanIm(i).im(:));
            
            % Remove NaNs 
            imTemp = imTemp(~isnan(imTemp));
            
            % Compute meadian and standard deviation of the mean pixel
            % brightness image
            imMedian=median(imTemp);
            imStd=std(imTemp);
            
            % Compute the brightness threshold for sun glint
            Glint(i).list=imStd*stdMag+imMedian;
            
            % Generate a logical for identifying pixels with brightness
            % above the glint threshold
            BinaryImage=meanIm(i).im>Glint(i).list;

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
            Glint_mask(i).list=statss(Indeks).PixelIdxList;                 % Old code: imStd*(stdMag+1)+imMedian;

        end
    end
end
