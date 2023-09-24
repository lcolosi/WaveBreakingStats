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
    % EO files. Function prints a warning if they are not consistent. 
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
        Indices(i)=str2num(arrTemp(k(2)+1:k(3)-1));
    
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
        
        % Save start and end indices for the ith track. 
        tracks_Im(i).Indices=[startI endI];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check for consistency between the raw imagery and EO trajectory file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: Check if the number of tracks and the number of time steps in each
    % track match between the raw imagery and EO file. 
    
    % Initialize boolean variable 
    True=1;
    
    %--- Check number of tracks ---% 
    if length(tracks)~=length(tracks_Im)
        True=0;
    end
    
    %--- Check number of images/EO file trajectory values ---% 
    if True==1
    
        % Loop through tracks
        for i=1:length(tracks)
    
            % Check if indices of the ith track from EO file is empty
            if ~isempty(tracks(i).Indices)
    
                % Compare the number of time steps between the EO file and
                % raw imagery
                if tracks_Im(i).Indices(2)-tracks_Im(i).Indices(1)~=tracks(i).Indices(2)-tracks(i).Indices(1)
                    True=0;
                end
            end
        end
    end
    
    % If tracks and images do not match, print warning statement
    if True==0
       fprintf('Number of tracks/images do not match!'); 
    end
