function tracks_Im = DefineTracksIm(D_Im,tracks)

    %%%%
    % tracks_Im = DefineTracksIm(D_Im,tracks)
    %
    % Function for identifying track numbers and determine indices of
    % each track using the raw imagery file names. A check for consistency
    % between the raw imagery and EO trajectory file is preformed. More
    % specifically, we check if the number of tracks and the number of
    % time steps (i.e., the number of images and the number for EO
    % trajectory values) in each track match between the raw imagery and
    % EO files. Also, it checks if the start and end file indicies match 
    % between raw imagery and EO files. Function throws an error if they
    % are not consistent. 
    %
    %   Parameters
    %   ----------
    %   D_Im     : Filenames of the non-georeferenced video images.
    %   tracks   : Structure containing the start and end time indices of
    %              each full flight track (located in the indices field).
    %              Indices derived from EO file. 
    % 
    %   Returns
    %   -------
    %   track_Im : Structure containing the start and end time indices of
    %              each full flight track (located in the indices field). 
    %              Indices derived from raw imagery file names. 
    % 
    %%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Identify track numbers and determine indices of tracks 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: Assume that the images are defined in form of :
    % CoreView_TrackNumber...
    
    % Initialize zeros array for track number ID (length = number of files)
    Indices=zeros(length(D_Im(:,1)),1);
    
    % Loop through nongeoreferenced files
    for i=1:length(D_Im(:,1))
    
        % Obtain file name of ith image
        arrTemp=D_Im(i,1).name;
    
        % Find the indices of first three underscores in the file name 
        k=find('_'==arrTemp,3,'first');
        
        % Obtain the track number ID between the second and third underscores
        Indices(i)=str2num(arrTemp(k(2)+1:k(3)-1));                         %#ok
    
    end
    
    % Reset the track numbers so they increment starting from 1
    Indices=Indices+1-min(Indices);
    
    % Obtain the highest track number
    maxI=max(Indices);
    
    % Loop through track numbers
    for i=1:maxI
    
        % Index of the first occurence of the ith track
        startI=find(i==Indices,1,'first');
        
        % Index of the last occurence of the ith track 
        endI=find(i==Indices,1,'last');
        
        % Save start and end time indicies for ith track if track is
        % greater than one time step
        if startI ~= endI
            tracks_Im(i).Indices=[startI endI];                             %#ok
        else
            tracks_Im(i).Indices = [];                                      %#ok
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for consistency between the raw imagery and EO trajectory file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: Check if the number of tracks and the number of time steps in each
    % track match between the raw imagery and EO file. 
    
    %--- Check number of tracks ---% 
    if length(tracks)~=length(tracks_Im)
        error('Number of tracks do not match!');
    end
    
    %--- Check number of images/EO file trajectory values and time indices ---% 
    
    % Loop through tracks
    for i=1:length(tracks)

        % Check if indices of the ith track from EO file is empty
        if ~isempty(tracks(i).Indices)

            % Compare the number of time steps between the EO file and
            % raw imagery
            if tracks_Im(i).Indices(2)-tracks_Im(i).Indices(1)~=tracks(i).Indices(2)-tracks(i).Indices(1)
                error('Number of images do not match!');
            end

            % Compare beginning and end time indices of the full track 
            % between images and EO files 
            if tracks_Im(i).Indices(1)~=tracks(i).Indices(1) || tracks_Im(i).Indices(2)~=tracks(i).Indices(2)
                error(['Start or end time indicies for Track ' num2str(i) ' between EO file and images do not match!']);
            end
        end
    end
end