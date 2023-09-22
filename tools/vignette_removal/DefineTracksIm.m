function tracks_Im = DefineTracksIm(D_Im,tracks)

    %%%%
    % tracks_Im = DefineTracksIm(D_Im,tracks)
    %
    % Function for 
    %
    %   Parameters
    %   ----------
    %   D_Im : Filenames of the non-georeferenced video images.
    %   tracks : Structure containing the start and end time indicies of
    %            each flight track.
    % 
    %   Returns
    %   -------
    %   track_Im : Structure with two fields: 
    %               (1) stable: An identifier for whether the Track is 
    %                           stable or unstable based on the inputted
    %                           roll, pitch, and heading criteria. Here, we
    %                           denote stability with the following
    %                           convention:
    %                                   0 -> Unstable
    %                                   1 -> Stable    
    %               (2) range: The start and end time indices of the stable
    %                          flight track.  
    %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine indices of tracks, and organize them. Assume that the images 
% are defined in form of : CoreView_TrackNumber.......

% 
Indices=zeros(length(D_Im(:,1)),1);

% Loop through nongeoreferenced files
for i=1:length(D_Im(:,1))

    % 
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
