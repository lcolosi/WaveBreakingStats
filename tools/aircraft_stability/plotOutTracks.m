function plotOutTracks(A,tracks,trackTag,dirstab,An,dirout)

    %%%%
    % plotOutTracks(A,tracks,trackTag,dirV,An,dirout)
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
    %   dirstab : Path to directory where figures will be saved.  
    %   An : Cropped GPSTime station cell array containing the time (UTC),
    %        track number, and image number.
    %   dirout : 
    % 
    %   Returns
    %   -------
    %   Three figures are saved to the dirstab directory: 
    %       (1) 
    %
    %   Question
    %   --------
    %   (1)   
    %
    %%%%

    % Set minimum extent of UTM Zone 10 coordinates
    minX=min(A.data(:,1));
    minY=min(A.data(:,2));
    
    hold on
    cm=jet(length(tracks));
    % Define initial time
    StartTime=datenum([2021,5,2,0,0,0]);
    load('D:\MASS\Processed\L_Computations\L_Computations\Test\TFOEx21_DEP01_BuoyGPS.mat');
    dirProcessed=[dirout 'Output\'];
    D_Im2=dir([dirProcessed 'Track_' num2str(1) '\' '*.tif']);
    proj=geotiffinfo([dirProcessed 'Track_' num2str(1) '\' D_Im2(11).name]);
    
    lonMax=-10^9;
    latMax=-10^9;
    lonMin=10^9;
    latMin=10^9;
    for i=1:length(tracks)
        indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);
        [latTemp,lonTemp] = projinv(proj,A.data(indicesPlot,1),A.data(indicesPlot,2));
        lonMax=max(lonMax,nanmax(lonTemp(:)));
        latMax=max(latMax,nanmax(latTemp(:)));
        lonMin=min(lonMin,nanmin(lonTemp(:)));
        latMin=min(latMin,nanmin(latTemp(:)));
    end
    lonMax=lonMax+0.1;
    latMax=latMax+0.1;
    lonMin=lonMin-0.1;
    latMin=latMin-0.1;
    % assignin('base','Temp1',lonMax);
    % assignin('base','Temp2',lonMin);
    for i=1:length(tracks)
        if (~isempty(tracks(i).Indices))%(trackTag(i).stable==1) & 
        figure('WindowState','maximized');
        
    
        vreme=str2num(cell2mat(An(tracks(i).Indices,1)));
        
        tracksTime=StartTime+(vreme(1)+vreme(2))/2/86400;
        lonn=interp1(BuoyTimeD1,-BuoyLON_D1,tracksTime);
        latt=interp1(BuoyTimeD1,BuoyLAT_D1,tracksTime);
        scatter(lonn,latt,'markerfacecolor',cm(i,:));
        if ~isempty(tracks(i).Indices)
            
            indicesPlot=tracks(i).Indices(1):1:tracks(i).Indices(2);
             [lat,lon] = projinv(proj,A.data(indicesPlot,1),A.data(indicesPlot,2));
            plot(lon,lat,'color',cm(i,:),'LineWidth',1.5);
            text(max(lon(1)-0.001),max(lat(1)-0.001),num2str(i),'color',cm(i,:),'FontName','times','fontsize',32)
    
    %         plot(A.data(indicesPlot,1)-minX,A.data(indicesPlot,2)-minY,'color',cm(i,:),'LineWidth',1.5);
    %         text(max(A.data(indicesPlot(end),1)-minX),max(A.data(indicesPlot(end),2)-minY),num2str(i),'color',cm(i,:),'FontName','times','fontsize',32)
        end
        end
        set(gca,'fontname','times','FontSize',32,'tickdir','both')
        xlabel('x (m)','fontname','times','FontSize',32,'Interpreter','latex')
        ylabel('y (m)','fontname','times','FontSize',32,'Interpreter','latex')
        title(['Track' num2str(i)],'fontname','times','FontSize',32,'Interpreter','latex')
        box on
        grid on
        xlim([lonMin lonMax])
        ylim([latMin latMax])
        saveas(gcf,[dirV 'Tracks_' num2str(i) '.png']);
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
        
        saveas(gcf,[dirV 'Track_Heading_' num2str(i) '.png']);
        close all
        end
end