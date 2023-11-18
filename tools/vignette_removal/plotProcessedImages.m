function plotProcessedImages(dirVn,dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,Glint,meanOriginal,outliers)

    %%%%
    % plotProcessedImages(dirVn,dirRaw,D_Im,tracks_Im,trackTag,tracks,meanIm,stdIm,Glint)
    %
    % Function for plotting the mean and standard Deviation Pixel 
    % Brightness image for each track and for plotting a side-by-side
    % comparison of the raw and processed (equilized/vignetting-removed)
    % images. 
    % 
    %
    %   Parameters
    %   ---------- 
    %   dirVn     : Path to directoryu for saving figures displaying
    %               vignette removal. 
    %   dirRaw    : Path to raw nongeoreferenced images. 
    %   D_Im      : Filenames of the raw nongeoreferenced video images.
    %   tracks_Im : Structure containing the start and end time indices of
    %               each full flight track (located in the indices field). 
    %               Indices derived from raw imagery file names. 
    %   trackTag  : Structure with three fields: 
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
    %   tracks   : Structure containing the start and end time indices of
    %              each full flight track (located in the indices field).
    %              Indices derived from EO file.
    %   meanIm   : Mean brightness of each pixel over the time period
    %              of the stable flight track (temporal mean). Mean
    %              brightness is smoothed using a 2D moving average.
    %   stdIm    : Sample standard deviation of brightness for each
    %               pixel over the time period of the stable flight 
    %               (temporal std). Standard deviation of brightness
    %               is smoothed using a 2D moving average.
    %   Glint    : Glint threshold brightness. 
    %   meanOriginal : Mean brightness over all pixels for each image 
    %                  (spatial average). Here, the mean of the all pixels
    %                  below the brightness threshold set by B_threshold of
    %                  each image is computed. This is computed to check
    %                  for significant variations in image lighting. 
    %   outliers : Statistics and theshold criteria for determining images
    %              with significant background brightness variations from 
    %              the temporal median bmean image brightness. 
    %
    %   Returns
    %   -------
    %   For each track, two sets of figures are returned: (1) mean and std
    %   of brightness for each pixel and (2) comparison between raw and
    %   processed images. 
    % 
    %%%%

    % Plotting parameters
    fontsize = 14; 
    font = 'times';  

    % Loop through tracks 
    for i=1 %1:length(tracks)
        
        % Check if track is stable and if full track has start and end time
        % indices
        if (trackTag(i).stable==1) && (~isempty(tracks(i).Indices))
            
            % Create a subdirectory for ith trajectory plots
            trackDir=[dirVn 'Track_' num2str(i) '\'];
            if ~exist(trackDir, 'dir'), mkdir(trackDir); end
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Plot the mean and std of brightness w.r.t. time for each pixel 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

            % Generate figure
            figure('Name', ['Mean and standard Deviation Pixel Brightness  - Flight Track' num2str(i)]);
            set(gcf,'color',[1,1,1])
            set(gcf,'Position',[100,100,1500,500])

            %--------- Subplot 1 ---------%
            subplot(1,2,1);
                
                % Plot mean pixel brightness 
                imagesc(meanIm(i).im);

                % set figure attributes
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                xlabel('x (pixels)','fontname',font,'FontSize',fontsize)
                ylabel('y (pixels)','fontname',font,'FontSize',fontsize)
                title(['Mean brightness - track ' num2str(i)],'fontname',font,'FontSize',fontsize)
                box on
                daspect([1 1 1])

                % Set colorbar attributes
                colormap(bone);
                yy=colorbar();
                title(yy,'Brightness','fontname',font,'FontSize',fontsize)

            %--------- Subplot 2 ---------%    
            subplot(1,2,2);

                % Plot standard deviation
                imagesc(stdIm(i).im);

                % set figure attributes
                set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                xlabel('x (pixels)','fontname',font,'FontSize',fontsize)
                ylabel('y (pixels)','fontname',font,'FontSize',fontsize)
                title(['Standard deviation of brightness - track ' num2str(i)],'fontname',font,'FontSize',fontsize)
                box on
                daspect([1 1 1])

                % Set colorbar attributes
                colormap(bone);
                yy=colorbar();
                title(yy,'Brightness','fontname',font,'FontSize',fontsize)
            
            % Save figure
            print(gcf,'-dpng', [trackDir 'Brightness_statistics-track' num2str(i) '.png'], '-r200');
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Plot the mean image brightness time series with outlier stats 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

            % Generate figure
            figure('Name', ['Mean Brightness per image  - Flight Track' num2str(i)]);
            set(gcf,'color',[1,1,1])
            set(gcf,'Position',[100,100,1200,500])
                
            % Plot mean image brightness (spatial average)
            plot(meanOriginal(i).nr, 'o-k', 'LineWidth',1);

            % Plot outlier statistics
            hold on 
                yline([outliers(i).nr.L outliers(i).nr.U outliers(i).nr.C],":",["Lower Threshold","Upper Threshold","Median Brightness"])
            hold off

            % set figure attributes
            set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
            xlabel('Number of frames','fontname',font,'FontSize',fontsize)
            ylabel('$\langle I(x,y,t) \rangle_{x,y}$','fontname',font,'FontSize',fontsize)
            box on
            %daspect([1 1 1])
        
            % Save figure
            print(gcf,'-dpng', [trackDir 'Mean_Image_Brightness_statistics-track' num2str(i) '.png'], '-r200');
    
            % Set beginning and end indices for stable flight period  
            beginSP = trackTag(i).range(1);
            endSP = trackTag(i).range(2);
            
            % Select 10 images from the track for plotting
            imNum=floor((0:9)*(endSP-beginSP)/10)+beginSP;

            % Loop through the ten images 
            for j=1:10
                
                % Load jth image
                a=imread([dirRaw D_Im(imNum(j)).name]);

                % Convert RGB image to grayscale and pixel values to double
                a=double(im2gray(a));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Crop high glint region of image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Compute the median of the mean and standard deviation images 
                meanTemp=median(meanIm(i).im(:),'omitnan');
                stdTemp=median(stdIm(i).im(:),'omitnan');
                
                % Remove mean pixel brightness and add spatial median of 
                % mean pixel brightness to equalize image 
                a2=a-meanIm(i).im+meanTemp;

                % Remove vignetting
                a2=mean(a2(:),'omitnan')+(a2-mean(a2(:),'omitnan'))./stdIm(i).im*stdTemp;
                
                % Obtain glint magnitude threshold
                glintMag=Glint(i).list;
                
                % Generate a logical array for identifying pixels with brightness
                % above the glint threshold
                BinaryImage = meanIm(i).im > glintMag;
                
                % Calculate the following properties for each region
                % with brightness greater than the glint threshold: 
                %   (1) Area: Number of pixels in the region (returns a
                %             scalar for each region identified).
                %   (2) PixelIdxList : Linear indices of the pixels in the
                %                     region (returns a vector of p 
                %                     elements long for each region where
                %                     p is the total number of pixels in
                %                     the region)
                % The statss variable is a Nx1 structure with fields Area
                % and PixelIdxList where N is the number of regions
                % identified.
                statss = regionprops(BinaryImage, 'Area','PixelIdxList');
                
                % Find the region with the largest area 
                [~,Indeks] = max( [statss.Area] );
                
                % Obtain the pixel indices of this region; this constitutes
                % the high glint mask
                GlintTemp=statss(Indeks).PixelIdxList;
                
                % Apply glint mask to image
                a2(GlintTemp)=nan;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Plot raw image and high-glint cropped image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
                
                % Generate figure
                figure('Name', ['Image Processing  - Flight Track' num2str(i)]);
                set(gcf,'color',[1,1,1])
                set(gcf,'Position',[100,100,1300,500])

                %--------- Subplot 1 ---------%
                subplot(1,2,1);
                    
                    % Plot raw image 
                    imagesc(a);

                    % set figure attributes
                    set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                    xlabel('x (pixels)','fontname',font,'FontSize',fontsize)
                    ylabel('y (pixels)','fontname',font,'FontSize',fontsize)
                    title(['Raw image - Image ' num2str(i)],'fontname',font,'FontSize',fontsize)
                    box on
                    daspect([1 1 1])

                    % Set colorbar attributes
                    colormap(bone);
                    yy=colorbar();
                    caxis([500 10000])
                    title(yy,'Brightness','fontname',font,'FontSize',fontsize)
                    
                %--------- Subplot 2 ---------%
                subplot(1,2,2);

                    % Plot high-glint cropped image
                    imagesc(a2);

                    % Set figure attributes
                    set(gca,'fontname',font,'FontSize',fontsize,'tickdir','both')
                    xlabel('x (pixels)','fontname',font,'FontSize',fontsize)
                    ylabel('y (pixels)','fontname',font,'FontSize',fontsize)
                    title(['High-glint cropped image - Image ' num2str(i)],'fontname',font,'FontSize',fontsize)
                    box on
                    daspect([1 1 1])

                    % Set colorbar attributes
                    colormap(bone);
                    yy=colorbar();
                    title(yy,'Brightness','fontname',font,'FontSize',fontsize)
                    caxis([500 10000])
                    drawnow

                % Save figure
                print(gcf,'-dpng', [trackDir 'Image_processing_compare-Image_M2' num2str(j) '.png'], '-r200');
        
            end
        end
        
    end
end