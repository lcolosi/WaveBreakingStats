function tracks = DefineTracks(An)

    %%%%
    % tracks = DefineTracks(An)
    %
    % Function for cutting out parts of the header that were not removed 
    % when importing data (look for '(sec)' in the string in the first 
    % column of A.textdata). The data field in the A structure is fine 
    % (A.data contains the position and attitude data). Only the textdata 
    % is cropped (A.textdata contains the GPSTime station). 
    %
    %   Parameters
    %   ----------
    %   A : Uncropped GPSTime station cell array.
    % 
    %   Returns
    %   -------
    %   A : Cropped GPSTime station cell array.            
    %
    %%%%

% Determine indices of tracks, and organize them. Assume that the files
% are defined in form of : ...._TrackNumber_ImageNumber .

Indices=zeros(length(An(:,1)),1);
for i=1:length(An(:,1))
    k=find('_'==An{i,2},2,'last');
    Indices(i)=str2num(An{i,2}(k(1)+1:k(2)-1));
end
Indices=Indices+1-min(Indices);
maxI=max(Indices);
for i=1:maxI
    startI=find(i==Indices,1,'first');
    endI=find(i==Indices,1,'last');
    tracks(i).Indices=[startI endI];
end