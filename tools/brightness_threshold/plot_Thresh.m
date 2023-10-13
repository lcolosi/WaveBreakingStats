function plot_Thresh(Threshold,N,dirProcessed,dirV,D_Im,tracks_Im,trackTag,tracks,arrayHist,Glint,Glint_mask,meanImm,GlobalOrLocal)
% Function for plotting images of mean brightness and std of brightness for each track,
% and for plotting out comparison of raw and images with removed vignetting
% Input - 

for i=11:12%1:length(tracks)
    if trackTag(i).stable==1
        trackDir=[dirV 'Track_' num2str(i) '\'];
        if ~exist(trackDir, 'dir'), mkdir(trackDir); end

        
        figure('WindowState','maximized');
        if GlobalOrLocal==1
            subplot(1,2,1);
            plot(N(i).hi(10:80),'color','black','linewidth',1.4);
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('Brightness','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('Counts','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Image histogram - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            grid on
            Array=gradient(gradient(log(movmean(N(i).hi,5))));
            subplot(1,2,2);
            plot(Array(10:80),'color','black','linewidth',1.4);
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('Brightness','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('Counts','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Second gradient of brightness - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            grid on
        elseif GlobalOrLocal==0
            subplot(1,2,1);
            imagesc(log(N(i).hi));
            colormap(jet(256))
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            ylabel('Image Nr.','fontname','times','FontSize',16,'Interpreter','latex')
            xlabel('Brightness','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Image histograms - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            grid on
            yy=colorbar();
            title(yy,'log(counts)','fontname','times','FontSize',16,'Interpreter','latex')
            
            subplot(1,2,2);
            plot(Threshold(i).th,'k')
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('Image Nr.','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('Brightness','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Breakness threshold - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            grid on
        end
        drawnow
        pathOut=trackDir;
        saveas(gcf,[pathOut 'MeanBrightnessStd.png'])
        %pathOut=trackDir;bdr_savefig(pathOut, 'MeanBrightnessStd', 'p', 600);
        close all
        % Select 10 images from the track for plotting
        beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
        endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);
        
        imNum=floor((0:9)*(endDif-beginDif)/10)+beginDif*0+1*1;
        
        D_Im=dir([dirProcessed 'Track_' num2str(i) '\' '*.tif']);
        
        %dirProcessed2=[dirProcessed 'Track_' num2str(i) '\'];
        dirProcessed2_mask=[dirProcessed 'TrackMask_' num2str(i) '\'];
        filenames_mask=dir([dirProcessed2_mask '*.tif']);
        for j=1:10
            
            % Load and plot raw image
            a=imread([dirProcessed 'Track_' num2str(i) '\' D_Im(imNum(j)).name]);
            
            a_mask=imread([dirProcessed2_mask D_Im(imNum(j)).name]);
            a_mask=double(im2gray(a_mask));
            %assignin('base','aaaaaaaaaaaa',meanImm);
            a_mask=a_mask*nanmax(meanImm(i).im(:))/256/256;
            BinaryImage=a_mask>Glint(i).list;
            statss = regionprops(BinaryImage, 'Area','PixelIdxList');
            [biggest_area,Indeks] = max( [statss.Area] );
            GlintTemp=statss(Indeks).PixelIdxList;
            a(GlintTemp)=nan;
            %assignin('base','aaaaaaaaaaaa',a);
            %assignin('base','aaaaaaaaaaaa2',a_mask);
            figure('WindowState','maximized');
            subplot(1,2,1);
            imagesc(a);
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Raw image - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            colormap(bone(256));
            daspect([1 1 1])
            initialI=nanmedian(a(:));
            if GlobalOrLocal==0
                caxis([5000 5001+Threshold(i).th(imNum(j))])
            else
                caxis([5000 5001+Threshold(i).th])
            end
            % Remove vignetting from raw image
            Both=1;
            
            
            %a2=RM_vig(a,meanIm(i).diff,stdIm(i).diff,meanIm(i).im,Both);
            subplot(1,2,2);
            
            if GlobalOrLocal==0
                imagesc(a>Threshold(i).th(imNum(j)));
            else
                imagesc(a>Threshold(i).th);
            end
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Processed image - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            colormap(bone(256));
            daspect([1 1 1])
            %caxis([3000 12000])
            drawnow
            pathOut=trackDir;
            saveas(gcf,[pathOut 'Comparison_' num2str(j) '.png'])
            %pathOut=trackDir;bdr_savefig(pathOut, ['Comparion_' num2str(j)], 'p', 600);
            close all
        end
        
    end
    
end