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
    %%%%

    % Minimum extent of UTM Zone 10 coordinates 
    minX=min(A.data(:,1));                                                  % Easting (Units: m)
    minY=min(A.data(:,2));                                                  % Northing (Units: m)

    % Plotting parameters
    cm=flipud(cbrewer2('RdYlBu',length(tracks)));                           % Colormap denoting the flight tracks.

    % Define initial time (needs to be an input argument)
    StartTime=datenum([2021,5,2,0,0,0]);

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
    %% Plot each individual flight track 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through tracks
    for i=1:length(tracks)

        % Check if time indicies for ith track are empty
        if ~isempty(tracks(i).Indices)

        % Generate figure
        figure('WindowState','maximized');
        
        % Obtain GPS time at beginning and end of ith full flight track    
        vreme=str2num(cell2mat(An(tracks(i).Indices,1)));
        
        % Convert to local time (incorrect calculation of local time; UTC time should be used)
        tracksTime=StartTime+(vreme(1)+vreme(2))/2/86400;
        
        % Obtain time indicies for ith flight track     
        indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);

        % Project UTM coorindates onto longitude and latitude grid (units: decimal degree)
        [lat,lon] = projinv(proj,A.data(indicesPlot,1),A.data(indicesPlot,2));
        
        % Plot ith track trajectory
        plot(lon,lat,'color',cm(i,:),'LineWidth',1.5);

        % Label ith flight track
        text(max(lon(1)-0.001),max(lat(1)-0.001),num2str(i),'color',cm(i,:),'FontName','times','fontsize',32)
        end

        % Set figure attributes
        set(gca,'fontname','times','FontSize',32,'tickdir','both')
        xlabel('x (m)','fontname','times','FontSize',32,'Interpreter','latex')
        ylabel('y (m)','fontname','times','FontSize',32,'Interpreter','latex')
        title(['Track' num2str(i)],'fontname','times','FontSize',32,'Interpreter','latex')
        box on
        grid on
        xlim([lonMin lonMax])
        ylim([latMin latMax])
        saveas(gcf,[dirTS 'Tracks_' num2str(i) '.png']);
        close all
    end
    
        
    
    
    for i=1:length(tracks)
        
        if ~isempty(tracks(i).Indices)
        if trackTag(i).stable==1
           answ='yes';
        else
            answ='no';
        end
        
        figure('WindowState','maximized');
        subplot(3,1,1)
        seriesPlot=A.data(tracks(i).Indices(1):tracks(i).Indices(2),6);
        dt=str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
        plot((1:length(seriesPlot))*dt,seriesPlot,'color','black','LineWidth',2)
        hold on
        % Indicates cropped out parts
        if trackTag(i).stable==1
            plot([trackTag(i).range(1)-tracks(i).Indices(1) trackTag(i).range(1)-tracks(i).Indices(1)]*dt,[-444 444],[trackTag(i).range(2)-tracks(i).Indices(1) trackTag(i).range(2)-tracks(i).Indices(1)]*dt,[-444 444],'color','red','LineWidth',2)
        end
        set(gca,'fontname','times','FontSize',32,'tickdir','both')
        %xlabel('Time (sec)','fontname','times','FontSize',32,'Interpreter','latex')
        ylabel('Heading ($^{\circ}$)','fontname','times','FontSize',32,'Interpreter','latex')
        title(['Angles - track ' num2str(i) ' Stable - ' answ],'fontname','times','FontSize',32,'Interpreter','latex')
        box on
        grid on
        ylim([mean(seriesPlot)-18 mean(seriesPlot)+18])
        xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
        
        subplot(3,1,2)
        seriesPlot=A.data(tracks(i).Indices(1):tracks(i).Indices(2),4);
        dt=str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
        plot((1:length(seriesPlot))*dt,seriesPlot,'color','black','LineWidth',2)
        hold on
        % Indicates cropped out parts
        if trackTag(i).stable==1
            plot([trackTag(i).range(1)-tracks(i).Indices(1) trackTag(i).range(1)-tracks(i).Indices(1)]*dt,[-444 444],[trackTag(i).range(2)-tracks(i).Indices(1) trackTag(i).range(2)-tracks(i).Indices(1)]*dt,[-444 444],'color','red','LineWidth',2)
        end
        set(gca,'fontname','times','FontSize',32,'tickdir','both')
        %xlabel('Time (sec)','fontname','times','FontSize',32,'Interpreter','latex')
        ylabel('Roll ($^{\circ}$)','fontname','times','FontSize',32,'Interpreter','latex')
        box on
        grid on
        ylim([mean(seriesPlot)-8 mean(seriesPlot)+8])
        xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
        
        
        subplot(3,1,3)
        seriesPlot=A.data(tracks(i).Indices(1):tracks(i).Indices(2),5);
        dt=str2num(An{tracks(i).Indices(1)+1,1})-str2num(An{tracks(i).Indices(1),1});
        plot((1:length(seriesPlot))*dt,seriesPlot,'color','black','LineWidth',2)
        hold on
        % Indicates cropped out parts
        if trackTag(i).stable==1
            plot([trackTag(i).range(1)-tracks(i).Indices(1) trackTag(i).range(1)-tracks(i).Indices(1)]*dt,[-444 444],[trackTag(i).range(2)-tracks(i).Indices(1) trackTag(i).range(2)-tracks(i).Indices(1)]*dt,[-444 444],'color','red','LineWidth',2)
        end
        set(gca,'fontname','times','FontSize',32,'tickdir','both')
        xlabel('Time (sec)','fontname','times','FontSize',32,'Interpreter','latex')
        ylabel('Pitch ($^{\circ}$)','fontname','times','FontSize',32,'Interpreter','latex')
        box on
        grid on
        ylim([mean(seriesPlot)-8 mean(seriesPlot)+8])
        xlim([0 tracks(i).Indices(2)-tracks(i).Indices(1)]*dt)
        
        saveas(gcf,[dirTS 'Track_Heading_' num2str(i) '.png']);
        close all
        end
    end

    %% Development Code


    load('D:\MASS\Processed\L_Computations\L_Computations\Test\TFOEx21_DEP01_BuoyGPS.mat');

    lonn=interp1(BuoyTimeD1,-BuoyLON_D1,tracksTime);
    latt=interp1(BuoyTimeD1,BuoyLAT_D1,tracksTime);
    scatter(lonn,latt,'markerfacecolor',cm(i,:));

    % assignin('base','Temp1',lonMax);
    % assignin('base','Temp2',lonMin);

