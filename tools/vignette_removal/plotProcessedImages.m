function plotProcessedImages(dirV,dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,Glint)
% Function for plotting images of mean brightness and std of brightness for each track,
% and for plotting out comparison of raw and images with removed vignetting
% Input - 

for i=11:12%1:length(tracks)
    if (trackTag(i).stable==1) & (~isempty(tracks(i).Indices))
        trackDir=[dirV 'Track_' num2str(i) '\'];
        if ~exist(trackDir, 'dir'), mkdir(trackDir); end
        
        
        figure('WindowState','maximized');
        subplot(1,2,1);
        imagesc(meanIm(i).im);
        set(gca,'fontname','times','FontSize',16,'tickdir','both')
        xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
        ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
        title(['Mean brightness - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
        box on
        daspect([1 1 1])
        colormap(bone(256));
        yy=colorbar();
        title(yy,'Brightness','fontname','times','FontSize',16)
        subplot(1,2,2);
        imagesc(stdIm(i).im);
        set(gca,'fontname','times','FontSize',16,'tickdir','both')
        xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
        ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
        title(['Mean brightness - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
        box on
        colormap(bone(256));
        yy=colorbar();
        title(yy,'Brightness','fontname','times','FontSize',16)
        daspect([1 1 1])
        drawnow
        
        pathOut=trackDir;
        saveas(gcf,[pathOut 'MeanBrightnessStd.png'])
        %bdr_savefig(pathOut, 'MeanBrightnessStd', 'p', 600);

        close all
        % Select 10 images from the track for plotting
        beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
        endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);
        
        imNum=floor((0:9)*(endDif-beginDif)/10)+beginDif;
        for j=1:10
            
            % Load and plot raw image
            a=imread([dirRaw D_Im(imNum(j)).name]);
            a=double(im2gray(a));
            
            
            
            figure('WindowState','maximized');
            subplot(1,2,1);
            imagesc(a);
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Raw image - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            colormap(jet(256));
            yy=colorbar();
            title(yy,'Brightness','fontname','times','FontSize',16)
            daspect([1 1 1])
            % Remove vignetting from raw image
            Both=1;
            
            meanTemp=median(meanIm(i).im(:));
            stdTemp=median(stdIm(i).im(:));
            
            a2=a-meanIm(i).im+meanTemp;
            a2=nanmean(a2(:))+(a2-nanmean(a2(:)))./stdIm(i).im*stdTemp;
            caxis([500 10000])
            %a2=RM_vig(a,meanIm(i).diff,stdIm(i).diff,meanIm(i).im,Both);
            subplot(1,2,2);
            
            glintMag=Glint(i).list;
            BinaryImage=meanIm(i).im>glintMag;
            statss = regionprops(BinaryImage, 'Area','PixelIdxList');
            [biggest_area,Indeks] = max( [statss.Area] );
            GlintTemp=statss(Indeks).PixelIdxList;
            a2(GlintTemp)=nan;
            
            imagesc(a2);
            set(gca,'fontname','times','FontSize',16,'tickdir','both')
            xlabel('x (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            ylabel('y (pixels)','fontname','times','FontSize',16,'Interpreter','latex')
            title(['Processed image - track ' num2str(i)],'fontname','times','FontSize',16,'Interpreter','latex')
            box on
            colormap(jet(256));
            daspect([1 1 1])
            caxis([500 10000])
            yy=colorbar();
            title(yy,'Brightness','fontname','times','FontSize',16)
            drawnow
            pathOut=trackDir;
            saveas(gcf,[pathOut 'Comparion_' num2str(j) '.png'])
            %bdr_savefig(pathOut, ['Comparion_' num2str(j)], 'p', 600);
            close all
        end
        
    end
    
end