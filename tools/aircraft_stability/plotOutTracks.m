function plotOutTracks(A,tracks,trackTag,dirTS,An,dirProc)

    %%%%
    % plotOutTracks(A,tracks,trackTag,dirTS,An,dirProc)
    %
    % Plotting function which preforms the following:  
    %   (1) Displays uncropped flight tracks, 
    %   (2) Indicating whether they are stable or not, 
    %   (3) Indicates which portion of the track was cropped (at the begining and
    %       at the end).
    %
    %   Parameters
    %   ----------
    %   A : Aircraft trajectory and atitude data array.
    %   tracks : Structure containing the start and end time indicies of
    %            each uncropped flight track.
    %   trackTag : Identifier for whether the Track is stable or unstable 
    %              based on the inputted roll, pitch, and heading criteria.
    %              Here, we denote stability with the following convention:
    %                               0 -> Unstable
    %                               1 -> Stable 
    %   dirTS : Path to directory where figures will be saved.  
    %   An : Cropped GPSTime station cell array containing the time (UTC),
    %        track number, and image number.
    %   dirProc : Path to Trimble processed (georeferenced) video images 
    %             for a given flight
    % 
    %   Returns
    %   -------
    %   Three figures are saved to the dirTS directory: 
    %       (1)  
    % 
    %   Notes
    %   -----
    %   This functions requires the mapping toolbox. Functions projinv and
    %   geotiffinto are apart of this toolbox. 
    % 
    %   Questions
    %   ---------
    %   Ask Nick about how to convert to GPS time and why there are weird
    %   features in the gps time variable. The method for deriving the time
    %   variable assumes that the data is collected at a regular interval,
    %   but there are definitely gaps in the time series which can cause
    %   the data to appear more spiky than it actually is. 
    %
    %%%%
    
    % Plotting parameters
    fontsize = 14; 
    font = 'times';
    cm=flipud(cbrewer2('RdYlBu',length(tracks)));                           % Colormap denoting the flight tracks. 
    gray = [0.5,0.5,0.5];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Obtain UTM coordinate info for projection to lat/lon grid 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set path to processed visible images files orgainzed by track
    dirProcessed=[dirProc 'Output\'];

    % Obtain filenames of of geotiff files fromtrack 1
    D_Im2=dir([dirProcessed 'Track_' num2str(1) '\' '*.tif']);

    % Grab image properties for one of the geotiff files
    proj=geotiffinfo([dirProcessed 'Track_' num2str(1) '\' D_Im2(1).name]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Find max and min lat and lon of all tracks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through tracks
    for i=1:length(tracks)

        % Obtain time indicies for ith flight track
        indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);

        % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
        [latTemp,lonTemp]=projinv(proj,A.data(indicesPlot,1),A.data(indicesPlot,2));

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
        [lat,lon] = projinv(proj,A.data(indicesPlot,1),A.data(indicesPlot,2));
        
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
    print(gcf,'-dpng', [dirTS 'flight_trjectory.png'], '-r300');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot aircraft stability for each track
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through tracks 
    for i=1:length(tracks)
        
        % Check if time indicies for ith track are empty 
        if ~isempty(tracks(i).Indices)

            % Check if cropped flight track is stable
            if trackTag(i).stable==1
               answ='Yes';
            else
                answ='No';
            end

            % Obtain time indicies for ith flight track     
            indicesFT=tracks(i).Indices(1):1:tracks(i).Indices(2);
            indicesCT = trackTag(i).range(1):trackTag(i).range(2); 

            % Grab roll, pitch, and heading time series data along entire
            % flight track
            h_ts = A.data(indicesFT,6); 
            p_ts = A.data(indicesFT,5);
            r_ts = A.data(indicesFT,4);

            % Compute a time vector for seconds since the begin of the
            % flight track
            dt = str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
            time_s = (1:length(h_ts))*dt; time_sc = (1:length(h_ts))*dt;

            % Generate figure
            figure('Name', ['Aircraft stability - Flight Track' num2str(i)]);
            set(gcf,'color',[1,1,1])
            set(gcf,'Position',[100,100,900,1300])

            %--------- Subplot 1 ---------%
            subplot(5,1,[1 2])

                % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
                [lat_ft,lon_ft] = projinv(proj,A.data(indicesFT,1),A.data(indicesFT,2));
                [lat_ct,lon_ct] = projinv(proj,A.data(indicesCT,1),A.data(indicesCT,2));
        
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
                plot(time_s,h_ts,'color','black','LineWidth',2)
                
                % Plot boundaries of stable flight
                hold on
                if trackTag(i).stable==1
                    % Compute the start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)
                end
    
                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                ylabel('Heading ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                % ylim([mean(h_ts)-8 mean(h_ts)+8])
                box on
                grid on
            
            %--------- Subplot 3 ---------%    
            subplot(5,1,4)
            
                % Plot from full flight trajectory
                plot(time_s,p_ts,'color','black','LineWidth',2)
                
                % Plot boundaries of stable flight
                hold on
                if trackTag(i).stable==1
                    % Compute the start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)
                end

                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                ylabel('Pitch ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                % ylim([mean(p_ts)-8 mean(p_ts)+8])
                box on
                grid on
            
            %--------- Subplot 4 ---------%
            subplot(5,1,5)
            
                % Plot from full flight trajectory
                plot(time_s,r_ts,'color','black','LineWidth',2)
                
                % Plot boundaries of stable flight
                hold on
                if trackTag(i).stable==1
                    % Compute the start and end times for the stable flight period
                    to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
                    tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
                    xline(to,'color','red','LineWidth',2)
                    xline(tf,'color','red','LineWidth',2)
                end

                % Set figure attributes 
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                xlabel('Time (sec)','fontname',font,'FontSize',fontsize)
                ylabel('Roll ($^{\circ}$)','fontname',font,'FontSize',fontsize)
                xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
                % ylim([mean(r_ts)-8 mean(r_ts)+8])
                box on
                grid on

            % Save figure
            print(gcf,'-dpng', [dirTS 'Stability_analysis-track' num2str(i) '.png'], '-r200'); 
        end
    end

    %% Development Code
    
%     load('D:\MASS\Processed\L_Computations\L_Computations\Test\TFOEx21_DEP01_BuoyGPS.mat');
% 
%     % Obtain GPS time at beginning and end of ith full flight track    
%     vreme=str2num(cell2mat(An(tracks(i).Indices,1)));
%     
%     % Convert to local time (incorrect calculation of local time; UTC time should be used)
%     tracksTime=StartTime+(vreme(1)+vreme(2))/2/86400;
% 
%     dt=str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
% 
%     lonn=interp1(BuoyTimeD1,-BuoyLON_D1,tracksTime);
%     latt=interp1(BuoyTimeD1,BuoyLAT_D1,tracksTime);
%     scatter(lonn,latt,'markerfacecolor',cm(i,:));
% 
%     % assignin('base','Temp1',lonMax);
%     % assignin('base','Temp2',lonMin);
% 
%     % Minimum extent of UTM Zone 10 coordinates 
%     minX=min(A.data(:,1));                                                  % Easting (Units: m)
%     minY=min(A.data(:,2));                                                  % Northing (Units: m)

% % Compute the start and end times for the stable flight period
% to = (trackTag(i).range(1)-tracks(i).Indices(1))*dt;
% tf = (trackTag(i).range(2)-tracks(i).Indices(1))*dt;
% 
% % Obtain the time series for the stable flight track 
% time_sc = to:dt:tf;