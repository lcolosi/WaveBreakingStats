function Threshold=denoiseThreshold(Threshold,tracks,trackTag)

for i=11:12%1:length(tracks)
    if (trackTag(i).stable==1) & (~isempty(tracks(i).Indices))
        [temp,indicesT]=rmoutliers(Threshold(i).th);
        temp=Threshold(i).th;
        temp(indicesT)=nan;
        while sum(double(isnan(temp)))>0
            temp=fillmissing(temp,'nearest');
        end
        Threshold(i).th=temp;
    end
    
end