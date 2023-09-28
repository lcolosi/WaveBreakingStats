function [Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag)

    %%%%
    % [Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag)
    %
    % Function for  
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
    %   Glint      : 
    %   Glint_mask :
    %  
    % 
    %%%%

    % Loop through tracks
    for i=2 %1:length(tracks)
        
        % Check if track is stable 
        if trackTag(i).stable==1
            
            % Sort mean pixel brightness in ascending numerical order
            imTemp=sort(meanIm(i).im(:));
            
            % 
            %imTemp=imTemp(1:floor(end/1));
            
            % Remove NaNs 
            imTemp = imTemp(~isnan(imTemp));
            
            % Compute meadian and standard 
            imMedian=median(imTemp);
            imStd=std(imTemp);
            
            % Find biggest connected area corresponding to glint using the
            % selected criteria.
            Glint(i).list=imStd*stdMag+imMedian;
            
            % Same, but for surrounding area
            Glint_mask(i).list=imStd*(stdMag+1)+imMedian;
        end
    end
end

%% Developmental code 
    % Find biggest connected area corresponding to glint using the
    % selected criteria.
    %         BinaryImage=meanIm(i).im>imMedian+imStd*stdMag;
    %         statss = regionprops(BinaryImage, 'Area','PixelIdxList');
    %         [biggest_area,Indeks] = max( [statss.Area] );
    %         Glint(i).list=statss(Indeks).PixelIdxList;