function WindDir=determine_wind(Speeds,WindTime,StartTime,tracks,An)

%ImageTime=StartTime+mean(str2num(cell2mat(An(tracks(1),1)))+str2num(cell2mat(An(tracks(2),1))))/86400/2;
ImageTime1=StartTime+0+(str2num(cell2mat(An(tracks(1),1))))/86400;
ImageTime2=StartTime+0+(str2num(cell2mat(An(tracks(2),1))))/86400;

Speeds(1)=10;
r = find(isnan(Speeds));
while sum(isnan(Speeds))>0
    
    Speeds(r) = Speeds(r-1);
end

tempN=interp1(WindTime,Speeds,ImageTime1:1/100000:ImageTime2);
%tempN=nanmean(tempN);
%assignin('base','aaaaaaaaaaaa',tempN);
% Note the vy is indexed with a negative sign!
WindDir.x=nanmean(cosd(tempN-180-90));
WindDir.y=nanmean(sind(tempN-180-90));
