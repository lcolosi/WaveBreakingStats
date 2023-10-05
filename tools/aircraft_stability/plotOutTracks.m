function plotOutTracks(A,tracks,trackTag,dirTS,An,sc,utc_time,StartDate,utmZone)

    %%%%
    % plotOutTracks(A,tracks,trackTag,dirTS,An,sc,utc_time,StartDate,utmZone)
    %
    % Plotting function which generates the following figures:  
    %   (1) The trajectory of all flight tracks on a lat-lon grid in movie
    %       and regular figure formats. 
    %   (3) The time intervals in UTC time of each track and its stable
    %       flight period. 
    %   (2) The individual flight trajectories with the time series of 
    %       roll, pitch, and heading.  
    % In addition, the titles of the third set of figures indicate whether
    % the track is stable or not.
    %
    %   Parameters
    %   ----------
    %   A        : Aircraft trajectory and atitude data array.
    %   tracks   : Structure containing the start and end time indicies of
    %              each uncropped flight track.
    %   trackTag : Structure with three fields: 
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
    %   dirTS    : Path to directory where figures will be saved.  
    %   An       : Cropped GPSTime station cell array containing the time 
    %              (UTC), track number, and image number.
    %   sc       : Stability criteria matrix which contains the criteria for the
    %              standard deviation for the roll, pitch, and heading (sigRoll,
    %              sigPitch, sigHeading respectively) as well as the number of
    %              standard deviations (Nstd) for the roll, pitch, or heading
    %              that constitute an abrupt change in attitude of the plane. sc 
    %              has the following array structure  
    %                     sc = [sigRoll,sigPitch,sigHeading,Nstd];
    %   utc_time  : UTC time in a datenum array.    
    %   StartDate : Start date string of flight in UTC in 'yyyymmdd' format. 
    %   utmZone   : UTM zone for experimental site as a string. 
    % 
    %   Returns
    %   -------
    %   All figures are saved to the dirTS directory.   
    %
    %%%%
    
    % Plotting parameters
    fontsize = 14; 
    font = 'times';
    cm=flipud(cbrewer2('RdYlBu',length(tracks)));                           % Colormap denoting the flight tracks. 
    gray = [0.5,0.5,0.5];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Find max and min lat and lon of all tracks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through tracks
    for i=1:length(tracks)

        % Check if time indicies for ith track are empty
        if ~isempty(tracks(i).Indices)

            % Obtain time indicies for ith flight track
            indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);
    
            % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
            [latTemp, lonTemp] = utm2deg(A.data(indicesPlot,1),A.data(indicesPlot,2), repmat(utmZone, length(A.data(indicesPlot,1)), 1));
    
            % Find longitude and latitude max and min for ith track and 
            % compare to previous flight tracks
            if i == 1
                lonMax=max(lonTemp,[],'all','omitnan'); latMax=max(latTemp,[],'all','omitnan');
                lonMin=min(lonTemp,[],'all','omitnan'); latMin=min(latTemp,[],'all','omitnan');
            else 
                lonMax=max([lonMax,max(lonTemp,[],'all','omitnan')]); latMax=max([latMax,max(latTemp,[],'all','omitnan')]);
                lonMin=min([lonMin,min(lonTemp,[],'all','omitnan')]); latMin=min([latMin,min(latTemp,[],'all','omitnan')]);
            end
        end
    end

    % Extend lon and lat maximum and minimum by 0.1 of a decimal degree
    lonMax=lonMax+0.1; latMax=latMax+0.1;
    lonMin=lonMin-0.1; latMin=latMin-0.1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot each individual flight track in a movie
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create object to write video files
    writerObj = VideoWriter([dirTS 'flight_trjectory.mp4'], 'MPEG-4');
    writerObj.FrameRate = 1;                                      
    open(writerObj);

    % Generate figure
    figure('Name', 'Flight Track and aircraft stability');
    set(gcf,'color',[1,1,1])
    set(gcf,'Position',[100,100,1300,900])

    % Loop through tracks
    for i=1:length(tracks)

        % Check if time indicies for ith track are empty
        if ~isempty(tracks(i).Indices)
        
            % Obtain time indicies for ith flight track     
            indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);
    
            % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
            [lat, lon] = utm2deg(A.data(indicesPlot,1),A.data(indicesPlot,2), repmat(utmZone, length(A.data(indicesPlot,1)), 1));
            
            % Plot ith track trajectory
            if i == 1
                geoscatter(lat,lon,30,cm(i,:),'o','filled');
            else 
                hold on 
                geoscatter(lat,lon,30,cm(i,:),'o','filled');
            end
    
            % Label ith flight track
            text(max(lat(1)+0.01),max(lon(1)-0.001),num2str(i),'color',cm(i,:),'FontName',font,'fontsize',fontsize)
            
            % Set figure attributes
            geobasemap satellite
            geolimits([latMin latMax],[lonMin lonMax])
            gx = gca; 
            gx.TickLabelFormat = '-dd';
            gx.LongitudeLabel.String = 'Longitude'; 
            gx.LatitudeLabel.String = 'Latitude';
            set(gx,'fontname',font,'FontSize',fontsize,'tickdir','both')
    
            % Set frame and write to video file
            frame = getframe(gcf);
            writeVideo(writerObj,frame)
        end
    end
    
    % Save Movie
    close(writerObj); 

    % Save figure 
    print(gcf,'-dpng', [dirTS 'flight_trajectory.png'], '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot time intervals associated with each track and its stable flight period 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate figure
    figure('Name', 'Aircraft stability Time Intervals','WindowState','maximized');
    set(gcf,'color',[1,1,1])

    % Loop through tracks
    for i=1:length(tracks)
        
        % Check if time indicies for ith track are empty 
        if ~isempty(tracks(i).Indices)
            
            % Obtain time indicies for ith flight track     
            indicesFT = tracks(i).Indices(1):1:tracks(i).Indices(2);
            indicesCT = trackTag(i).range(1):trackTag(i).range(2);

            % Plot ith track time interval for full and stable flight
            % periods
            hold on
                plot(indicesFT,utc_time(indicesFT),'.','Color',gray,'MarkerSize',8) 
                plot(indicesCT,utc_time(indicesCT),'.','Color',cm(i,:),'MarkerSize',8)

            % Label ith flight track
            text(max(indicesFT-1100),max(utc_time(indicesFT)),num2str(i),'color','k','FontName',font,'fontsize',fontsize)
        end
    end

    % Set figure attributes 
    set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
    xlabel('Number of measurements','fontname',font,'FontSize',fontsize)
    ylabel(['UTC time from ' StartDate(end-3:end-2) '/' StartDate(end-1:end) '/' StartDate(1:4) '(hrs)'],'fontname',font,'FontSize',fontsize)
    datetick('y','HH','keeplimits')
    box on
    grid on

    % Save figure 
    print(gcf,'-dpng', [dirTS 'Stable_flight_periods.png'], '-r300');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot aircraft stability for each track
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through tracks 
    for i= 1:length(tracks)
        
        % Check if time indicies for ith track are empty 
        if ~isempty(tracks(i).Indices)

            % Check if cropped flight track is stable
            if trackTag(i).stable==1
               answ='Yes';
            else
                answ='No';
            end

            % Obtain time indicies for ith flight track     
            indicesFT = tracks(i).Indices(1):1:tracks(i).Indices(2);
            indicesCT = trackTag(i).range(1):trackTag(i).range(2); 

            % Grab roll, pitch, and heading time series data along entire
            % flight track
            h_ts = A.data(indicesFT,6); 
            p_ts = A.data(indicesFT,5);
            r_ts = A.data(indicesFT,4);

            % Compute a time vector for seconds since the begin of the
            % flight track
            dt = str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
            time_s = (1:length(h_ts))*dt; 

            % Generate figure
            figure('Name', ['Aircraft stability - Flight Track' num2str(i)]);
            set(gcf,'color',[1,1,1])
            set(gcf,'Position',[100,100,900,1300])

            %--------- Subplot 1 ---------%
            subplot(5,1,[1 2])

                % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
                [lat_ft,lon_ft] = utm2deg(A.data(indicesFT,1),A.data(indicesFT,2), repmat(utmZone, length(A.data(indicesFT,1)), 1));
                [lat_ct,lon_ct] = utm2deg(A.data(indicesCT,1),A.data(indicesCT,2), repmat(utmZone, length(A.data(indicesCT,1)), 1));
        
                % Get max and mins and add 0.05 decimal degrees
                lonmaxTemp = max(lon_ft,[],'all','omitnan')+0.01; latmaxTemp = max(lat_ft,[],'all','omitnan')+0.01;
                lonminTemp = min(lon_ft,[],'all','omitnan')-0.01; latminTemp = min(lat_ft,[],'all','omitnan')-0.01;

                % Plot ith track trajectory
                geoscatter(lat_ft,lon_ft,30,gray,'o','filled');
                hold on 
                geoscatter(lat_ct,lon_ct,30,'red','o','filled');

                % Set figure attributes
                title(['Track ' num2str(i) ' | Stable - ' answ],'fontname',font,'FontSize',fontsize)
                geobasemap satellite
                geolimits([latminTemp, latmaxTemp], [lonminTemp, lonmaxTemp])
                gx = gca; 
                gx.TickLabelFormat = '-dd';
                gx.LongitudeLabel.String = 'Longitude'; 
                gx.LatitudeLabel.String = 'Latitude';
                set(gx,'fontname',font,'FontSize',fontsize,'tickdir','both')

            %--------- Subplot 2 ---------%
            subplot(5,1,3)
                
                % Plot from full flight trajectory
                plot(time_s,h_ts,'.','color','black','MarkerSize',5)
                
                % Plot boundaries of stable flight and stability criteria
                hold on
                if trackTag(i).stable==1
                    % Display start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)
                    
                    % Display stability criteria
                    max_std = sc(4)*trackTag(i).stats(3,2)+trackTag(i).stats(3,1);
                    min_std = -sc(4)*trackTag(i).stats(3,2)+trackTag(i).stats(3,1);
                    yline([max_std min_std],'--b',{'$t\sigma_h$', '$-t\sigma_h$'})

                    % Display initial stable range (determined from maxPer
                    % and shift)
                    idx_o = tracks(i).Indices(1)+trackTag(i).stats(4,1);    
                    idx_f = tracks(i).Indices(2)-trackTag(i).stats(4,2);
                    time_si = time_s(indicesFT >= idx_o & indicesFT <= idx_f);
                    fx = [time_si(1), time_si(1), time_si(end), time_si(end)];
                    fy = [min_std-4, max_std+4, max_std+4, min_std-4];
                    fill(fx,fy,'b','FaceAlpha',0.1,'EdgeColor','none')
                end
    
                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                ylabel('Heading ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                ylim([min_std-4 max_std+4])                                
                box on
                grid on
            
            %--------- Subplot 3 ---------%    
            subplot(5,1,4)
            
                % Plot from full flight trajectory
                plot(time_s,p_ts,'.','color','black','MarkerSize',5)
                
                % Plot boundaries of stable flight
                hold on
                if trackTag(i).stable==1
                    % Compute the start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)

                    % Display stability criteria
                    max_std = sc(4)*trackTag(i).stats(2,2)+trackTag(i).stats(2,1);
                    min_std = -sc(4)*trackTag(i).stats(2,2)+trackTag(i).stats(2,1);
                    yline([max_std min_std],'--b',{'$t\sigma_p$', '$-t\sigma_p$'})

                    % Display initial stable range (determined from maxPer
                    % and shift)
                    idx_o = tracks(i).Indices(1)+trackTag(i).stats(4,1); 
                    idx_f = tracks(i).Indices(2)-trackTag(i).stats(4,2);
                    time_si = time_s(indicesFT >= idx_o & indicesFT <= idx_f);
                    fx = [time_si(1), time_si(1), time_si(end), time_si(end)];
                    fy = [min_std-4, max_std+4, max_std+4, min_std-4];
                    fill(fx,fy,'b','FaceAlpha',0.1,'EdgeColor','none')
                end

                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                ylabel('Pitch ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                ylim([min_std-4 max_std+4]) 
                box on
                grid on
            
            %--------- Subplot 4 ---------%
            subplot(5,1,5)
            
                % Plot from full flight trajectory
                plot(time_s,r_ts,'.','color','black','MarkerSize',5)
                
                % Plot boundaries of stable flight
                hold on
                if trackTag(i).stable==1
                    % Compute the start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)

                    % Display stability criteria
                    max_std = sc(4)*trackTag(i).stats(1,2)+trackTag(i).stats(1,1);
                    min_std = -sc(4)*trackTag(i).stats(1,2)+trackTag(i).stats(1,1);
                    yline([max_std min_std],'--b',{'$t\sigma_r$', '$-t\sigma_r$'})

                    % Display initial stable range (determined from maxPer
                    % and shift)
                    idx_o = tracks(i).Indices(1)+trackTag(i).stats(4,1); 
                    idx_f = tracks(i).Indices(2)-trackTag(i).stats(4,2);
                    time_si = time_s(indicesFT >= idx_o & indicesFT <= idx_f);
                    fx = [time_si(1), time_si(1), time_si(end), time_si(end)];
                    fy = [min_std-4, max_std+4, max_std+4, min_std-4];
                    fill(fx,fy,'b','FaceAlpha',0.1,'EdgeColor','none')
                end

                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                xlabel('Time (sec)','fontname',font,'FontSize',fontsize)
                ylabel('Roll ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                ylim([min_std-4 max_std+4]) 
                box on
                grid on

            % Save figure
            print(gcf,'-dpng', [dirTS 'Stability_analysis-track' num2str(i) '.png'], '-r200'); 
        end
    end
end
