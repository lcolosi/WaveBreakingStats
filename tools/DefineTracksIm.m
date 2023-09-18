function tracks_Im = DefineTracksIm(D_Im,tracks)

% Determine indices of tracks, and organize them. Assume that the images 
% are defined in form of : CoreView_TrackNumber.......

Indices=zeros(length(D_Im(:,1)),1);
for i=1:length(D_Im(:,1))
    arrTemp=D_Im(i,1).name;
    k=find('_'==arrTemp,3,'first');
    Indices(i)=str2num(arrTemp(k(2)+1:k(3)-1));
end
Indices=Indices+1-min(Indices);
maxI=max(Indices);
for i=1:maxI
    startI=find(i==Indices,1,'first');
    endI=find(i==Indices,1,'last');
    tracks_Im(i).Indices=[startI endI];
end

% Check if the number of tracks, and number of images in each track match
% between the EO file and the raw images.
True=1;
if length(tracks)~=length(tracks_Im)
    True=0;
end
if True==1
    for i=1:length(tracks)
        if ~isempty(tracks(i).Indices)
            if tracks_Im(i).Indices(2)-tracks_Im(i).Indices(1)~=tracks(i).Indices(2)-tracks(i).Indices(1)
                True=0;
            end
        end
    end
end

if True==0
   fprintf('Number of tracks/images not matching!'); 
end
