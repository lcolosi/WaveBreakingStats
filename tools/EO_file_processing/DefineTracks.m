function tracks = DefineTracks(An)

    %%%%
    % tracks = DefineTracks(An)
    %
    % Function for determining the start and end time indices of tracks 
    % (i.e., flight track), and relabel flight tracks starting at 1. 
    % Here, we assumes that the files are defined in form of: 
    %               MASS_VIDEO_TrackNumber_ImageNumber .
    %
    %   Parameters
    %   ----------
    %   An     : Cropped GPSTime station cell array containing the time (UTC),
    %            track number, and image number. 
    % 
    %   Returns
    %   -------
    %   tracks : Structure containing the start and end time indicies of
    %            each flight track. 
    %
    %%%%

    % Set indices array with a length equal to the number of time steps in
    % the time series. 
    Indices=zeros(length(An(:,1)),1);

    % Loop through time 
    for i=1:length(An(:,1))

        % Obtain the indices of the underscores infront of the track 
        % number and image number 
        k=find('_'==An{i,2},2,'last');

        % Obtain the track number for the ith time step 
        Indices(i)=str2num(An{i,2}(k(1)+1:k(2)-1));
    end

    % Make the track numbers increment starting from 1
    Indices=Indices+1-min(Indices);
    
    % Find the highest track number 
    maxI=max(Indices);
    
    % Loop through track numbers
    for i=1:maxI

        % Obtain the start and end time indicies of the ith flight line 
        startI=find(i==Indices,1,'first');
        endI=find(i==Indices,1,'last');

        % Save start and end time indicies in tracks structure if track is
        % greater than one time step
        if startI ~= endI
            tracks(i).Indices=[startI endI];
        else
            tracks(i).Indices = [];
        end
    end

end