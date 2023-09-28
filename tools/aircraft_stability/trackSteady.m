function trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,Shift,t,tCheck,An)

    %%%%
    % trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,,shift,t,tcheckAn)
    %
    % Function for definning whether tracks is steady by analyzing roll,
    % pitch, and heading of aircraft. In order to avoid possible change of
    % direction at the end/start of the track, the beginning and end of 
    % the track are removed from consideration.
    % 
    % Additionally, the function determines (based on the relative change
    % of angles) which portion of the beginning and/or end of the flight
    % track needs to be cropped.
    %
    %   Parameters
    %   ----------
    %   A          : Aircraft trajectory and atitude data array.
    %   tracks     : Structure containing the start and end time indicies of
    %                each flight track.
    %   maxPer     : Maximum percent of flight track to be removed at the begin
    %                and end of the track. For example, maxPer = [25,25] means
    %                that at a minimium, 50 percent (from 25% - 75%) of the
    %                flight track will be used in analysis with the 0-25%
    %                and 75%-100% removed from the flight track.
    %   sigRoll    : Maximum allowed roll standard deviations in the stable
    %                segment. 
    %   sigPitch   : Maximum allowed pitch standard deviations in the stable
    %                segment.
    %   sigHeading : Maximum allowed heading standard deviations in the 
    %                stable segment.
    %   Shift      : Array of size Nx2 containing N potenital shifts in 
    %                beginning and end percentages of the flight track in case
    %                the flight track cropped by maxPer is unstable. Column 1
    %                contains shifts in beginning flight track index and Column 2
    %                contains shifts in end flight track index. 
    %   t          : Number of standard deviations of either roll, pitch, 
    %                or heading that constitute an abrupt change in 
    %                attitude of the plane.
    %   tcheck     : Time interval between the jth roll/pitch/heading 
    %                observation to check if there is an abrupt change in 
    %                attitude shortly after the jth observation. 
    %                Units: seconds
    %   An         : Cropped GPSTime station cell array containing the time 
    %                (UTC), track number, and image number.
    % 
    %   Returns
    %   -------
    %   trackTag   : Structure with three fields: 
    %               (1) stable: An identifier for whether the Track is 
    %                           stable or unstable based on the inputted
    %                           roll, pitch, and heading criteria. Here, we
    %                           denote stability with the following
    %                           convention:
    %                                   0 -> Unstable
    %                                   1 -> Stable    
    %               (2) range: The start and end time indices of the stable
    %                          flight track. 
    %               (3) stats: An array containing the statistics from the
    %                          initial stability analysis using the maxPer 
    %                          and shift parameters. The array has the
    %                          following structure: 
    %                                  [mean_r, std_r;... 
    %                                   mean_p, std_p;...
    %                                   mean_h, std_h; 
    %                                   rmN(1), rmN(2)]
    %                           where mean and standard deviations are
    %                           computed over the interval: 
    %                                   tracks(i).Indices(1)+rmN(1) to 
    %                                   tracks(i).Indices(2)-rmN(2)
    %                           
    %%%%

    % Loop through tracks
    for i=1:length(tracks)

        % Check if the start and end time indices for the ith track are empty 
        if ~isempty(tracks(i).Indices)
            
            % Set stablity of ith track to zero (denoting unstable) for the
            % case where the three time interval shifts do not find a stable
            % flight track
            trackTag(i).stable=0;

            % Set maximum percent of flight track to be removed at the
            % begin and end of the track
            maxPerTemp=maxPer;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Check for stable flight period using maxPer and shift
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Loop through three time interval shifting cases if flight
            % track with beginning and end cropped is unstable 
            for jjj=1:size(Shift,1)

                % Set maximum number of images to be removed at begining/end of flight track
                rmN=floor(maxPerTemp/100*(tracks(i).Indices(2)-tracks(i).Indices(1)+1));
                
                % Set the roll, pitch, and heading from the EO file for the selected track
                roll=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),4);
                pitch=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),5);
                heading=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),6);
                
                % Set variables for computing directional stats 
                dt = NaN; task = false;                                     % Parameters for computing standard error of the mean

                % Compute standard deviation of roll, pitch, and heading
                mean_r = mean(roll); std_r = std(roll);                     
                mean_p = mean(pitch); std_p = std(pitch);
                [mean_h,std_h, ~] = direction_stats(heading, dt, task);     % Forces heading to vary from 0 to 360 degrees

                % Check stability criteria for each attitude angle
                if std_r>sigRoll || std_p>sigPitch || std_h>sigHeading
                    % Flight Track is unstable
                    trackTag(i).stable=0;
                else
                    % Flight track is stable
                    trackTag(i).stable=1;
                end

                % Case 1: Flight track is stable -> end loop 
                if trackTag(i).stable==1
                    break

                % Case 2:Flight track is not stable -> try shifting initial
                % time period by shift.
                else
                    maxPerTemp=maxPer+Shift(jjj,:);
                end
                
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Final determine time indices for stable flight track 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % After verifying flight stability by the standard deviation of
            % the entire cropped flight track, we can look to extend the 
            % stable flight until we hit an abrupt change in angle. Here, 
            % we will increment toward beginning and end of the full flight
            % track and check each time step for the following: 
            %
            % (1) The jth roll/pitch/heading observation deviates from the
            %     its mean by more than t times its standard deviation 
            %     plus 1 (for roll and pitch) or 2 (for heading).    
            %
            % (2) The jth roll/pitch/heading observation 
            %     7 seconds away (when the jth observation is within 7
            %     seconds, use the time index at the beginning/end of the
            %     track) deviates from its mean by more than t times its
            %     standard deviation plus 1 (for roll and pitch) or 2 
            %     (for heading).
            % 
            % Once, we hit an abrupt change in angle or are approaching a
            % time period with an abrupt direction change, the loop breaks
            % and the index is saved.  
            % 
            % Note that we use the mean and standard deviation of the
            % roll/pitch/heading from the stable time period found above by
            % applying maxPer and shift parameters.  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Compute the time interval between measurements (units: seconds)
            dt = str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
            
            %--------- Roll ---------%
            % Loop through time indices for beginning of flight track in decending order  
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                
                % Set minimum time index based on roll
                minR=j;

                % Check if jth observation is an abrupt change in roll
                if abs(A.data(j,4)-mean_r)>t*std_r && abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),4)-mean_r)>t*std_r
                    break
                end
            end

            % Loop through time indicies for end if flight track in ascending order 
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                
                % Set maximum time index based on roll
                maxR=j;

                % Check if jth observation is an abrupt change in roll
                if abs(A.data(j,4)-mean_r)>t*std_r && abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),4)-mean_r)>t*std_r
                    break
                end
            end
            
            %--------- Pitch ---------%
            % Loop through time indices for beginning of flight track in decending order  
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                
                % Set minimum time index based on pitch
                minP=j;
                
                % Check if jth observation is an abrupt change in pitch
                if abs(A.data(j,5)-mean_p)>t*std_p && abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),5)-mean_p)>t*std_p
                    break
                end
            end

            % Loop through time indicies for end if flight track in ascending order 
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                
                % Set maximum time index based on pitch
                maxP=j;

                % Check if jth observation is an abrupt change in pitch
                if abs(A.data(j,5)-mean_p)>t*std_p && abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),5)-mean_p)>t*std_p
                    break
                end
            end

            %--------- Heading ---------%
            % Loop through time indices for beginning of flight track in decending order  
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                
                % Set minimum time index based on heading
                minH=j;

                % Check if jth observation is an abrupt change in heading
                if abs(A.data(j,6)-mean_h)>t*std_h && abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),6)-mean_h)>t*std_h
                    break
                end
            end

            % Loop through time indicies for end if flight track in ascending order 
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                
                % Set maximum time index based on heading
                maxH=j;

                % Check if jth observation is an abrupt change in heading
                if abs(A.data(j,6)-mean_h)>t*std_h && abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),6)-mean_h)>t*std_h
                    break
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Save flight stability results in trackTag 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Flight track is unstable -> range set to NaN 
            if trackTag(i).stable==0
                trackTag(i).range=nan;
                trackTag(i).stats=nan;

            % Flight track is stable -> set lower limit to maximum time 
            % index of minimum roll, pitch, and heading time indices (vice
            % versa for upper limit) and save statistics
            else
                trackTag(i).range=[max([minR,minP,minH]) min([maxR,maxP,maxH])];
                trackTag(i).stats=[mean_r, std_r;... 
                                   mean_p, std_p;...
                                   mean_h, std_h; 
                                   rmN(1), rmN(2)];
            end
        
        % If the start and end time indices for the ith track are empty, label the track as unstable 
        else
            trackTag(i).stable=0;
        end
    end
end