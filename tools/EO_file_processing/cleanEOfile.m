function cleanEOfile(EOpath,EOdir,tracks,D)

    %%%%
    % cleanEOfile(EOpath,EOdir,tracks,D)
    %
    % Function for creating EO files for each track and saves them to the
    % directory of the parent EO file. 
    %
    %   Parameters
    %   ---------- 
    %   EOpath : Path to EO file in EOdir 
    %   EOdir  : Path to directory where EO text are located. 
    %   tracks : Structure containing the start and end time indicies of
    %            each flight track. 
    %   D      : Filenames of all EO files.
    % 
    %   Returns
    %   -------
    %   Creates EO files for each track and saves them to the directory of
    %   the parent EO file. The parent EO file is the original EO 
    %   containing data from the entire flight.  
    %%%%

    % Generate waitbar
    f = waitbar(0,'Please wait...');

    % Loop through tracks
    for i=1:length(tracks)

        % Check if the start and end time indices for the ith track are empty 
        if (~isempty(tracks(i).Indices))

            % Update waitbar
            waitbar(i/length(tracks),f,...
                {['Running cleanEOfile.m : On track ' num2str(i) ' of ' num2str(length(tracks))]; [num2str(round((i/length(tracks))*100)) '$\%$ complete...']})

            % Generate EO filename for the ith track and append it EO 
            % directory path 
            ind = strfind(D(1).name,'_');
            tempStr=[EOdir D(1).name(1:ind(end)) 'track_' num2str(i) '_' D(1).name(1+ind(end):end)];
            
            % Open parent and new EO files 
            fidN = fopen(tempStr,'w');
            fid=fopen(EOpath);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Copy header over to new EO file for the ith track
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Read first line of parent EO and print it in new EO file
            tline = fgetl(fid);
            fprintf(fidN,'%s\n',tline);

            % Set counter
            Counter=1;

            % Loop through lines in text file and check if line contains '(sec)'
            % string (marks the end of the header)
            while ~contains(tline,'(sec)')
                % Print ith line in parent EO to new EO file
                tline = fgetl(fid);
                fprintf(fidN,'%s\n',tline);

                % Reset counter
                Counter=Counter+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Copy data from the ith track to new EO file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Read parent EO file and search for lines not containing a new
            % line (i.e. lines containing data)
            file=fileread(EOpath);
            lines = regexp(file, '\n', 'split');
            
            % Loop through time steps for ith track and print data in new
            % EO file
            for j=tracks(i).Indices(1):tracks(i).Indices(2)
                fprintf(fidN,'%s\n',lines{1,j+Counter});
            end
            
            % Close EO files
            fclose(fid);
            fclose(fidN);
        end
    end
    
    % Close waitbar
    close(f)

end