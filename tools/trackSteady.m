function trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,An)

    %%%%
    % trackTag = trackSteady(A,tracks,maxPer,sigRoll,sigPitch,sigHeading,An)
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
    %   A : Aircraft trajectory and atitude data array.
    %   tracks : Structure containing the start and end time indicies of
    %            each flight track.
    %   maxPer : Step size in search.
    %   sigRoll : Maximum allowed roll standard deviations in the stable
    %             segment. 
    %   sigPitch : Maximum allowed pitch standard deviations in the stable
    %              segment.
    %   sigHeading : Maximum allowed heading standard deviations in the 
    %                stable segment.
    %   An : Cropped GPSTime station cell array containing the time (UTC),
    %        track number, and image number.
    % 
    %   Returns
    %   -------
    %   trackTag :  
    %
    %%%%


    % Loop through tracks
    for i=1:length(tracks)

        % Check if the start and end time indices for the ith track are empty 
        if ~isempty(tracks(i).Indices)
            
            % Set stablity of ith track to zero (denoting unstable)
            trackTag(i).stable=0;

            % Set maximum 
            maxPerTemp=maxPer;

            % 
            Shift=[-0.8*maxPer(1) 0.8*maxPer(2);0.8*maxPer(1) -0.8*maxPer(2);0 0];
            
            % Check if there is a region where the flight is steady
            for jjj=1:3
                % number of images to be removed at begining/end
                rmN=floor(maxPerTemp/100*(tracks(i).Indices(2)-tracks(i).Indices(1)+1));
                
                % Load roll, pitch, and heading from the EO file for the selected track
                roll=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),4);
                pitch=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),5);
                heading=A.data(tracks(i).Indices(1)+rmN(1):tracks(i).Indices(2)-rmN(2),6);
                
                if mean(heading)>270
                    heading(heading<90)=heading(heading<90)+360;
                elseif mean(heading)<90
                    heading(heading>270)=heading(heading>270)-360;
                end
                
                if std(roll)>sigRoll | std(pitch)>sigPitch | std(heading)>sigHeading
                    trackTag(i).stable=0;
                else
                    trackTag(i).stable=1;
                end
                
                % If the track is not stable try shifting initial region of
                % interest by Shift.
                if trackTag(i).stable==1
                    break
                else
                    maxPerTemp=maxPer+Shift(jjj,:);
                end
                
            end
            
            dt=str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
            
            % determine range of indices. If the track is not stable return nan
            
            stdVar=4;
            % compute min and max indice based on roll. These are determined as a
            % point where the value varies from the mean by more than stdVar times
            % the standard deviation.
            
            % roll
            meanR=mean(roll);
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                minR=j;
                % In order to avoid brief changes in angles to cause a break, check
                % if the same is true for point on track 7 (tCheck) seconds apart
                tCheck=7;
                if abs(A.data(j,4)-meanR)>4*std(roll)+1 & abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),4)-meanR)>4*std(heading)+1
                    break
                end
            end
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                maxR=j;
                if abs(A.data(j,4)-meanR)>4*std(roll)+1 & abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),4)-meanR)>4*std(heading)+1
                    break
                end
            end
            
            %pitch
            meanP=mean(pitch);
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                minP=j;
                if abs(A.data(j,5)-meanP)>4*std(pitch)+1 & abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),5)-meanP)>4*std(heading)+1
                    break
                end
            end
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                maxP=j;
                if abs(A.data(j,5)-meanP)>4*std(pitch)+1 & abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),5)-meanP)>4*std(heading)+1
                    break
                end
            end
            %heading
            meanH=mean(heading);
            for j=tracks(i).Indices(1)+rmN(1)-1:-1:tracks(i).Indices(1)
                minH=j;
                if abs(A.data(j,6)-meanH)>4*std(heading)+2 & abs(A.data(max(tracks(i).Indices(1),j-floor(tCheck/dt)),6)-meanH)>4*std(heading)+2
                    break
                end
            end
            for j=tracks(i).Indices(2)-rmN(2)+1:1:tracks(i).Indices(2)
                maxH=j;
                if abs(A.data(j,6)-meanH)>4*std(heading)+2 & abs(A.data(min(tracks(i).Indices(2),j+floor(tCheck/dt)),6)-meanH)>4*std(heading)+2
                    break
                end
            end
            
            % define ranges
            if trackTag(i).stable==0
                trackTag(i).range=nan;
            else
                trackTag(i).range=[max([minR,minP,minH]) min([maxR,maxP,maxH])];
            end
        else
            trackTag(i).stable=0;
        end
    end
