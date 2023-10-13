function Threshold=denoiseThreshold(Threshold,tracks,trackTag)

    %%%%
    % Threshold=denoiseThreshold(Threshold,tracks,trackTag)
    %
    % Function for      
    %
    %   Parameters
    %   ---------- 
    %   Threshold : 
    %   tracks    : 
    %   trackTag  : 
    % 
    %   Returns
    %   ------- 
    %   Threshold : 
    % 
    %%%%

    % Loop through tracks
    for i= 1:length(tracks)


        if (trackTag(i).stable==1) && (~isempty(tracks(i).Indices))


            [~,indicesT]=rmoutliers(Threshold(i).th);


            temp=Threshold(i).th;
            temp(indicesT)=nan;
            
            
            while sum(double(isnan(temp)))>0
                temp=fillmissing(temp,'nearest');
            end
            
            
            Threshold(i).th=temp;
        end
    end
end