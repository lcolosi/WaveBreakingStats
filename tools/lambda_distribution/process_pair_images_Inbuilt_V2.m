function [Lambda,Lambda2,Area_Tot] = process_pair_images_Inbuilt_V2(A,B,C,D,PARAM_dt,PARAM_dx,alphaa,Threshold,Thresh_t,Thresh_x,WindDir,SpreadMax,iop,dirin,filenames2)
%% Setup
warning off
addpath('utilities')
addpath('optical_flow')

%Temp
% im1=A;
% im2=B;

% Set space and time resolution
dt = PARAM_dt;
dx = PARAM_dx;
%th = PARAM_threshold_contour; % threshold number for whitecap detection (grayscale)
th=2;
box_padding = 15;  %pixels;

% Optical flow settings
alpha = alphaa;
ratio = 0.85;
minWidth = 20;
nOuterFPIterations = 15;
nInnerFPIterations = 1;
nSORIterations = 40;

% get image area. Note that txe pixels near the edges (same as box_padding
% are discarded from the considerations.
ImTemp=A(:,:,1);
ImTemp(ImTemp==0)=nan;
ImTemp=conv2(ImTemp,ones(box_padding,box_padding)/box_padding/box_padding,'same');
Area_Tot=sum(A(~isnan(ImTemp)),'all');

ImTemp=A(:,:,1);
ImTemp2=A(:,:,1);
ImTemp2(ImTemp2==0)=nan;
% test removing area around edges
lengthW=floor(5/PARAM_dx)*2;
ImTemp2=conv2(ImTemp2,ones(lengthW,lengthW),'same');
ImTemp(isnan(ImTemp2))=0;


%ImTemp=A(:,:,1);
ImTemp(ImTemp<Threshold)=0;
% convert image to binary image
bw = im2bw(uint16(ImTemp),th./2^16);
%assignin('base','AA',bw);
% temp
%save original image
% Bb=B;

% Determine stats for all areas
[Bound,L] = bwboundaries(bw,'noholes');

stats = regionprops(L,'Area','Centroid','MajorAxisLength' ,'MinorAxisLength','orientation','perimeter','pixellist');
for i=1:length(stats)
    Area(i) = stats(i).Area;
end
siz=max(0.4/dx^2,7);
[selected_regions] = find(Area>siz);

% The optical flow algorithm
para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
I1 = A;
I2 = B;

%assignin('base','stats',stats)

%% Determine velocities
k=0;
% length(selected_regions)
for counter=1:length(selected_regions)
%     length(selected_regions)-counter
    % find box encompasssing region
    pixel_list = stats(selected_regions(counter)).PixelList;
    
    % Bounds of areas
    min_px_x = min(pixel_list(:,1));
    max_px_x = max(pixel_list(:,1));
    min_px_y = min(pixel_list(:,2));
    max_px_y = max(pixel_list(:,2));
    
    % Images for which the velocities are computed
    check1=1;



    Size=2^ceil(log2(max([1+max_px_x-min_px_x+2*box_padding 1+max_px_y-min_px_y+2*box_padding])));
    stepXmin=floor((Size-(1+max_px_x-min_px_x))/2);
    stepXmax=ceil((Size-(1+max_px_x-min_px_x))/2);
    stepYmin=floor((Size-(1+max_px_y-min_px_y))/2);
    stepYmax=ceil((Size-(1+max_px_y-min_px_y))/2);
    
    if min_px_x-stepXmin-5>0 && min_px_y-stepYmin-5>0 && max_px_y+stepYmax+5<length(A(:,1,1)) && max_px_x+stepXmax+5<length(A(1,:,1))
        subImageA = I1(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
        subImageB = I2(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
        subImageC = C(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
        subImageD = D(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
        if sum(double((subImageA(:))==0))>0 || sum(double((subImageB(:))==0))>0 || sum(double((subImageC(:))==0))>0 || sum(double((subImageD(:))==0))>0
            check1=0;
        end
    else
        check1=0;
    end
%     subImageA = I1(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
%     subImageB = I2(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
%     subImageC = C(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);
%     subImageD = D(min_px_y-stepYmin:max_px_y+stepYmax,min_px_x-stepXmin:max_px_x+stepXmax,:);

    
    
    if (check1==1)
        k=k+1;
        
%         if mod(k,100)==0
%            k 
%         end
        pixel_list_local(:,2) = pixel_list(:,2)-(min_px_y-stepYmin)+1;
        pixel_list_local(:,1) = pixel_list(:,1)-(min_px_x-stepXmin)+1;

        % contour
        pixel_list_contour = Bound{selected_regions(counter)};
        pixel_list_contour2(:,1) = pixel_list_contour(:,2);
        pixel_list_contour2(:,2) = pixel_list_contour(:,1);
        pixel_list_contour = pixel_list_contour2;
        clear pixel_list_contour2;
        
        pixel_list_local_contour(:,2) = pixel_list_contour(:,2)-(min_px_y-stepYmin)+1;
        pixel_list_local_contour(:,1) = pixel_list_contour(:,1)-(min_px_x-stepXmin)+1;
        Lambda(k).ContourList_px =  pixel_list_contour;
        
        subImageA=double(subImageA);
        subImageB=double(subImageB);
        subImageC=double(subImageC);
        subImageD=double(subImageD);
        
        % Optical flow coefficient (25% of mean brightness at the contour)
        im1=A(:,:,1);
        pix1=0;
        lengthPix=0;
        for i=1:length(Lambda(k).ContourList_px)
            ind1=Lambda(k).ContourList_px(i,2);
            ind2=Lambda(k).ContourList_px(i,1);
            pix1=pix1+double(im1(ind1,ind2));
            lengthPix=lengthPix+1;
        end
        alpha=pix1/lengthPix/16;
        % test alternative
        alpha=mean(subImageA(:)-5000)/32;
        alpha=(Threshold-5000)/24;
        meanBright=subImageA(pixel_list_local(:,2),pixel_list_local(:,1));
        alpha=(mean(meanBright(:))-5000)/24;
        
        % get velocity - Optical flow technique
        
        % Assign appropriate area size. This essentially determines the
        % appropriate initial downsampling. For larger breakers 6x6 pixel
        % area is represented with one pixel. This enables for easier
        % determination of large velocities, however downsampling small
        % breakers to that size often leads to errors.
        sizeArea=min(size(subImageA))-2*box_padding;
        if sizeArea<4
            sizeArea=3;
        elseif (4<sizeArea) && (sizeArea<6)
            sizeArea=4;
        elseif (6<sizeArea) && (sizeArea<8)
            sizeArea=5;   
        else
            sizeArea=6;
        end
        
        minWidth=floor(size(subImageA,2)/sizeArea);
        
        % Determine velocities. In the optical flow method, the velocities
        % are constrained with the alpha parameter (penalizes large
        % velocity gradients in the energy functional). As the same values
        % are not appropriate for breakers of different brightness, the
        % results are iteratively checked for large outliers. The criteria
        % are: velocities above 20 m/s, ratio of speeds of neighboring
        % pixels larger than 200%, and velocity gradients larger than 6 m/s
        % at neighboring pixels. 
        counter2=0;
        while counter2<1
            counter2=2;
            para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
            image1=subImageA(:,:,1);
            image2=subImageB(:,:,1);
            image3=subImageC(:,:,1);
            image4=subImageD(:,:,1);

            maxx1=max(image1(:));
            maxx2=max(image2(:));
            maxx3=max(image3(:));
            maxx4=max(image4(:));
            maxx=max([maxx1 maxx2 maxx3 maxx4]);
            image1=image1/maxx;
            image2=image2/maxx;
            image3=image3/maxx;
            image4=image4/maxx;
            
            
            
            %% Scale images
%             testIm1=linspace(0.9,1.1,length(image1(1,:)));
%             testIm3=linspace(1.1,0.9,length(image1(:,1)));
%             
%             
%             image1=image1.*testIm1.*testIm3';
%             image2=image2.*testIm1.*testIm3';
            
            %% First velocity
            ImgStack(:,:,1) = image1;% combine the images to a sequence
            ImgStack(:,:,2) = image2;
            if exist('opticFlow')
                reset(opticFlow); % clear previous optical flow
            end
            opticFlow = opticalFlowFarneback('numpyramidlevels',4,'pyramidscale',0.5,'filtersize',11,'NumIterations',15,'NeighborhoodSize',3);
            for tOptical = 1:2
                flow = estimateFlow(opticFlow,ImgStack(:,:,tOptical));
            end
            vx=flow.Vx;
            vy=flow.Vy;

            %% Second velocity
            ImgStack2(:,:,1) = image2;% combine the images to a sequence
            ImgStack2(:,:,2) = image3;
            if exist('opticFlow')
                reset(opticFlow); % clear previous optical flow
            end
            opticFlow = opticalFlowFarneback('numpyramidlevels',4,'pyramidscale',0.5,'filtersize',11,'NumIterations',15,'NeighborhoodSize',3);
            for tOptical = 1:2
                flow = estimateFlow(opticFlow,ImgStack2(:,:,tOptical));
            end
            vxC=flow.Vx;
            vyC=flow.Vy;

            %% Third velocity
            ImgStack3(:,:,1) = image3;% combine the images to a sequence
            ImgStack3(:,:,2) = image4;
            if exist('opticFlow')
                reset(opticFlow); % clear previous optical flow
            end
            opticFlow = opticalFlowFarneback('numpyramidlevels',4,'pyramidscale',0.5,'filtersize',11,'NumIterations',15,'NeighborhoodSize',3);
            for tOptical = 1:2
                flow = estimateFlow(opticFlow,ImgStack3(:,:,tOptical));
            end
            vxD=flow.Vx;
            vyD=flow.Vy;
            

            
            for i=1:length(pixel_list_local_contour(:,1))
                vx_region2(i) = vx(pixel_list_local_contour(i,2),pixel_list_local_contour(i,1));
                vy_region2(i) = vy(pixel_list_local_contour(i,2),pixel_list_local_contour(i,1));
            end
%             for i=1:length(pixel_list_local_contour(:,1))
%                 vxC_region2(i) = vxC(pixel_list_local_contour(i,2),pixel_list_local_contour(i,1));
%                 vyC_region2(i) = vyC(pixel_list_local_contour(i,2),pixel_list_local_contour(i,1));
%             end
            
            Lambda(k).VX_contour_m_s =  vx_region2*dx/dt;
            Lambda(k).VY_contour_m_s =  vy_region2*dx/dt;
            [GridX,GridY]=meshgrid(1:size(subImageA,2),1:size(subImageA,1));


            %% Position of pixels in subsequent images
            xB_temp=pixel_list_local_contour(:,1)+vx_region2';
            yB_temp=pixel_list_local_contour(:,2)+vy_region2';
            vx_temp = interp2(GridX,GridY,vxC,xB_temp,yB_temp);
            vy_temp = interp2(GridX,GridY,vyC,xB_temp,yB_temp);
            BoundaryB=interp2(GridX,GridY,subImageB(:,:,1),xB_temp,yB_temp);
            xC_temp=xB_temp+vx_temp;
            yC_temp=yB_temp+vy_temp;
            vxC_temp = interp2(GridX,GridY,vxD,xC_temp,yC_temp);
            vyC_temp = interp2(GridX,GridY,vyD,xC_temp,yC_temp);
            BoundaryC=interp2(GridX,GridY,subImageC(:,:,1),xC_temp,yC_temp);

            xD_temp=xC_temp+vxC_temp;
            yD_temp=yC_temp+vyC_temp;
            BoundaryD=interp2(GridX,GridY,subImageD(:,:,1),xD_temp,yD_temp);
%             for i=1:length(pixel_list_local_contour(:,1))
%                 xC_temp(i)=pixel_list_local_contour(i,1)+vx_region2(i);
%                 yC_temp(i)=pixel_list_local_contour(i,2)+vy_region2(i);
%                 vx_temp(i) = interp2(GridX,GridY,vxC,xC_temp(i),yC_temp(i));
%                 vy_temp(i) = interp2(GridX,GridY,vyC,xC_temp(i),yC_temp(i));
%                 
% %                 length1=length(subImageC(:,1,1));
% %                 length2=length(subImageC(1,:,1));
% %                 
% %                 cTemp1=subImageC(max(1,floor(yC_temp(i))),max(1,floor(xC_temp(i))),1);
% %                 cTemp2=subImageC(max(1,floor(yC_temp(i))),min(length2,ceil(xC_temp(i))),1);
% %                 cTemp3=subImageC(min(length1,ceil(yC_temp(i))),max(1,floor(xC_temp(i))),1);
% %                 cTemp4=subImageC(min(length1,ceil(yC_temp(i))),min(length2,ceil(xC_temp(i))),1);
% %                 BoundaryC(i)=max([cTemp1,cTemp2,cTemp3,cTemp4]);
%                 BoundaryC(i)=interp2(GridX,GridY,subImageC(:,:,1),xC_temp(i),yC_temp(i));
%             end
%             for i=1:length(pixel_list_local_contour(:,1))
%                 xD_temp(i)=xC_temp(i)+vx_temp(i);
%                 yD_temp(i)=yC_temp(i)+vy_temp(i);
%                 vxD_temp(i) = interp2(GridX,GridY,vxD,xD_temp(i),yD_temp(i));
%                 vyD_temp(i) = interp2(GridX,GridY,vyD,xD_temp(i),yD_temp(i));
%                 BoundaryD(i)=interp2(GridX,GridY,subImageD(:,:,1),xD_temp(i),yD_temp(i));
%             end
            
            
            
            Lambda(k).ContourLocal=pixel_list_local_contour;
            Lambda(k).xC=xC_temp;
            Lambda(k).yC=yC_temp;
            Lambda(k).vxB_temp=vx_temp;
            Lambda(k).vyB_temp=vy_temp;
            Lambda(k).xD=xD_temp;
            Lambda(k).yD=yD_temp;
            Lambda(k).vxC_temp=vxC_temp;
            Lambda(k).vyC_temp=vyC_temp;

            Lambda(k).xB=xB_temp;
            Lambda(k).yB=yB_temp;
            clear xC_temp yC_temp vx_temp vy_temp xB_temp yB_temp xD_temp yD_temp vxC_temp vyC_temp ImgStack ImgStack2 ImgStack3
            
         
            
%             maxDistance=0;
%             maxRatio=0;
% %             for qq=1:length(Lambda(k).ContourList_px(:,1));
% %                 [maxV,maxR]=findNeighbours(qq,Lambda(k).ContourList_px,Lambda(k).VX_contour_m_s,Lambda(k).VY_contour_m_s);
% %                 maxDistance=max(maxDistance,maxV);
% %                 maxRatio=max(maxRatio,maxR);
% %             end
% %             
% %             if (maxDistance>6) | (maxRatio>1) | (max(sqrt( Lambda(k).VX_contour_m_s.^2 + Lambda(k).VY_contour_m_s.^2 ))>20)
% %                 alpha=alpha*1.3;
% %                 counter2=counter2+1;
% % 
% %             else
% %                 % arbitrary large value that breaks the where loop
% %                 counter2=111;
% %             end
%             counter2=111;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Lambda(k).subImageA=subImageA;
%         Lambda(k).subImageB=subImageB;
%         Lambda(k).subImageC=subImageC;
%         Lambda(k).subImageD=subImageD;
        

%         if length(Lambda(k).xC)>30
%             [MaxV,IndeksMax]=sort(sqrt(Lambda(k).VX_contour_m_s.^2+Lambda(k).VY_contour_m_s.^2));
%             Test1=sum(sqrt((vx_region2(IndeksMax(end-50:end))-Lambda(k).vxC_temp(IndeksMax(end-50:end))).^2+(vy_region2(IndeksMax(end-50:end))-Lambda(k).vyC_temp(IndeksMax(end-50:end))).^2))./ ...
%             sum(sqrt((vx_region2(IndeksMax(end-50:end))).^2+(vy_region2(IndeksMax(end-50:end))).^2));
%             Test2=sum(sqrt((Lambda(k).vxD_temp(IndeksMax(end-50:end))-Lambda(k).vxC_temp(IndeksMax(end-50:end))).^2+(Lambda(k).vyD_temp(IndeksMax(end-50:end))-Lambda(k).vyC_temp(IndeksMax(end-50:end))).^2))./ ...
%             sum(sqrt((Lambda(k).vxC_temp(IndeksMax(end-50:end))).^2+(Lambda(k).vyC_temp(IndeksMax(end-50:end))).^2));
%         else
%             Test1=sum(sqrt((vx_region2-Lambda(k).vxC_temp).^2+(vy_region2-Lambda(k).vyC_temp).^2))./ ...
%             sum(sqrt((vx_region2).^2+(vy_region2).^2));
%             Test2=sum(sqrt((Lambda(k).vxD_temp-Lambda(k).vxC_temp).^2+(Lambda(k).vyD_temp-Lambda(k).vyC_temp).^2))./ ...
%             sum(sqrt((Lambda(k).vxC_temp).^2+(Lambda(k).vyC_temp).^2));
%         end

        %% Check consistency of velocities
        if length(Lambda(k).xC)>40
            [MaxV,IndeksMax]=sort(sqrt(Lambda(k).VX_contour_m_s.^2+Lambda(k).VY_contour_m_s.^2));
            Test1C=corr(vx_region2(IndeksMax(end-40:end))',Lambda(k).vxB_temp(IndeksMax(end-40:end)));
            Test2C=corr(vy_region2(IndeksMax(end-40:end))',Lambda(k).vyB_temp(IndeksMax(end-40:end)));

        else
            Test1C=corr(vx_region2',Lambda(k).vxB_temp);
            Test2C=corr(vy_region2',Lambda(k).vyB_temp);

        end

        if length(Lambda(k).xC)>40
            [MaxV,IndeksMax]=sort(sqrt(Lambda(k).VX_contour_m_s.^2+Lambda(k).VY_contour_m_s.^2));
            Test1D=corr(Lambda(k).vxC_temp(IndeksMax(end-40:end)),Lambda(k).vxC_temp(IndeksMax(end-40:end)));
            Test2D=corr(Lambda(k).vyC_temp(IndeksMax(end-40:end)),Lambda(k).vyC_temp(IndeksMax(end-40:end)));

        else
            Test1D=corr(Lambda(k).vxC_temp,Lambda(k).vxB_temp);
            Test2D=corr(Lambda(k).vyC_temp,Lambda(k).vyB_temp);

        end
        
        
        Lambda(k).testC=(Test1C+Test2C)/2;
        Lambda(k).testD=(Test1D+Test2D)/2;
        Lambda(k).test2C=sum(BoundaryC>Threshold-50)/length(BoundaryC);
        Lambda(k).test2D=sum(BoundaryD>Threshold-50)/length(BoundaryD);

        %% Number of pixels in subimageD above threshold

        xD=Lambda(k).xD;
        yD=Lambda(k).yD;
        Contour_at_D_x=[];
        Contour_at_D_y=[];
        for iT=1:length(Lambda(k).xD)
            if iT==1
                Contour_at_D_x(iT)=xD(1);
                Contour_at_D_y(iT)=yD(1);
                xD(iT)=[];
                yD(iT)=[];
            else
                [Md,Id]=min(sqrt((xD-Contour_at_D_x(iT-1)).^2+(yD-Contour_at_D_y(iT-1)).^2));
                Contour_at_D_x(iT)=xD(Id);
                Contour_at_D_y(iT)=yD(Id);
                xD(Id)=[];
                yD(Id)=[];
            end

        end


        Lambda(k).Contour_at_D_x = Contour_at_D_x ;
        Lambda(k).Contour_at_D_y = Contour_at_D_y;
        Lambda(k).Area_m2 =  stats(selected_regions(counter)).Area*dx*dx;
        Lambda(k).Area_D_m2=polyarea(Contour_at_D_x,Contour_at_D_y)*dx*dx;

        [xD_all,yD_all]=meshgrid(1:length(subImageD(1,:,1)),1:length(subImageD(:,1,1)));

        [IN]    = inpolygon(xD_all(:),yD_all(:),Contour_at_D_x,Contour_at_D_y);
        IN      = double(IN);
        subD=subImageD(:,:,1);
        subD=subD(:);
        subD=subD(IN==1);
        
        Lambda(k).Area_D_WC_m2=sum(double(subD>Threshold))*dx*dx;



        clear BoundaryC BoundaryD
        if 2>1%(Lambda(k).test>0.65)&&(Lambda(k).test2>0.3)
            %% Determining active breakers and assigning values 
            for i=1:length(pixel_list_local(:,1))
                vx_region(i) = vx(pixel_list_local(i,2),pixel_list_local(i,1));
                vy_region(i) = vy(pixel_list_local(i,2),pixel_list_local(i,1));
            end


            % Assigning variables
            Lambda(k).VX_m_s =  vx_region*dx/dt;
            Lambda(k).VY_m_s =  vy_region*dx/dt;
            
            Lambda(k).PixelList_px =  stats(selected_regions(counter)).PixelList;
            Lambda(k).Centroid_px =  stats(selected_regions(counter)).Centroid;

            % This variable (tttt)
            tttt=zeros(1,length(Lambda(k).VX_contour_m_s));


            ttt=[];
            ttt2=[];
            %% Space gradient
            tttt_X_grad=zeros(1,length(Lambda(k).VX_contour_m_s));
            X_grad=zeros(1,length(Lambda(k).VX_contour_m_s));
            [fx,fy]=gradient(double(A(:,:,1)));
            for i=1:length(Lambda(k).VX_contour_m_s)
                ind1=Lambda(k).ContourList_px(i,2);
                ind2=Lambda(k).ContourList_px(i,1);
                X_grad(i)=sqrt(fx(ind1,ind2)^2+fy(ind1,ind2)^2);
                if sqrt(fx(ind1,ind2)^2+fy(ind1,ind2)^2)>Thresh_x
                    tttt(i)=tttt(i)+1;
                    tttt_X_grad(i)=1;
                end
            end
            im2=B(:,:,1);
            
            %% time gradient
            tttt_T_grad=zeros(1,length(Lambda(k).VX_contour_m_s));
            T_grad=zeros(1,length(Lambda(k).VX_contour_m_s));
            for i=1:length(Lambda(k).VX_contour_m_s)
                ind1=Lambda(k).ContourList_px(i,2);
                ind2=Lambda(k).ContourList_px(i,1);
                pix1=double(im2(ind1,ind2));
                pix2=double(im1(ind1,ind2));
                T_grad(i)=abs(pix1-pix2);
                if abs(pix1-pix2)>Thresh_t
                    tttt(i)=tttt(i)+1;
                    tttt_T_grad(i)=1;
                end
            end

            tttt_DIR=zeros(1,length(Lambda(k).VX_contour_m_s));
            % Direction criteria 
            for i=1:length(Lambda(k).VX_contour_m_s)
                [theta,Rn] = cart2pol(Lambda(k).VX_contour_m_s(i),Lambda(k).VY_contour_m_s(i));

                u=[Lambda(k).VX_contour_m_s(i) Lambda(k).VY_contour_m_s(i)];

                % unit vector for wind
                v=[WindDir.x WindDir.y];
                CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
                ThetaInDegrees = real(acosd(CosTheta));

                if ThetaInDegrees<=SpreadMax%(theta<=-1.8326 | theta>=0.2618)
                    tttt(i)=tttt(i)+1;
                    tttt_DIR(i)=1;
                end
            end

            mean_VXn = nanmean(Lambda(k).VX_contour_m_s);
            mean_VYn = nanmean(Lambda(k).VY_contour_m_s);

            tttt_Nor=zeros(1,length(Lambda(k).VX_contour_m_s));
            % Determine normals
            curvNall=[];
            lenEl=[];
            for i=1:length(Lambda(k).VX_contour_m_s)
                if length(Lambda(k).VX_contour_m_s)>3
                    aa1=Lambda(k).ContourList_px(:,1)-Lambda(k).ContourList_px(i,1);
                    aa2=Lambda(k).ContourList_px(:,2)-Lambda(k).ContourList_px(i,2);
                    tre=sqrt(aa1.^2+aa2.^2);
                    [nNorm,In]=mink(tre,5);

                    treN=tre;
                    treN(treN==0)=22;
                    [nNormN,InN]=mink(treN,2);

                    lenEl=[lenEl sum(nNormN)/2];

                    if max(aa1(In))-min(aa1(In))>=max(aa2(In))-min(aa2(In))
                        fNorm=polyfit(Lambda(k).ContourList_px(In,1)-min(Lambda(k).ContourList_px(In,1)),Lambda(k).ContourList_px(In,2)-min(Lambda(k).ContourList_px(In,2)),1);
                        ug=fNorm(1)*(Lambda(k).ContourList_px(i,1)-min(Lambda(k).ContourList_px(In,1)))+fNorm(2);
                        %assignin('base','x',Lambda(k).ContourList_px(In,1)-min(Lambda(k).ContourList_px(In,1)));
                        %assignin('base','y',Lambda(k).ContourList_px(In,2)-min(Lambda(k).ContourList_px(In,2)));
    %                     fNorm=fit(Lambda(k).ContourList_px(In,1)-min(Lambda(k).ContourList_px(In,1)),Lambda(k).ContourList_px(In,2)-min(Lambda(k).ContourList_px(In,2)),'poly1');
    %                     ug=fNorm.p1*(Lambda(k).ContourList_px(i,1)-min(Lambda(k).ContourList_px(In,1)))+fNorm.p2;
    %                     ug=fNorm.p1;
                        %ug=fNorm.p1*(Lambda(k).ContourList_px(i,1)-min(Lambda(k).ContourList_px(In,1)))^2+fNorm.p2*(Lambda(k).ContourList_px(i,1)-min(Lambda(k).ContourList_px(In,1)))+fNorm.p3;
                        curvN=[ug -1];
                        curvN=curvN/max(abs(curvN))/2;
                        if inpolygon(Lambda(k).ContourList_px(i,1)+curvN(1),Lambda(k).ContourList_px(i,2)+curvN(2), Lambda(k).ContourList_px(:,1)', Lambda(k).ContourList_px(:,2)')
                            curvN=[-ug 1];
                            curvN=curvN/max(abs(curvN))/2;
                            if inpolygon(Lambda(k).ContourList_px(i,1)+curvN(1),Lambda(k).ContourList_px(i,2)+curvN(2), Lambda(k).ContourList_px(:,1)', Lambda(k).ContourList_px(:,2)')
                                curvN=[-1 -1];


                            end
                        end
                    else
                        fNorm=polyfit(Lambda(k).ContourList_px(In,2)-min(Lambda(k).ContourList_px(In,2)),Lambda(k).ContourList_px(In,1)-min(Lambda(k).ContourList_px(In,1)),1);
                        ug=fNorm(1)*(Lambda(k).ContourList_px(i,2)-min(Lambda(k).ContourList_px(In,2)))+fNorm(2);

    %                     fNorm=fit(Lambda(k).ContourList_px(In,2)-min(Lambda(k).ContourList_px(In,2)),Lambda(k).ContourList_px(In,1)-min(Lambda(k).ContourList_px(In,1)),'poly1');
    %                     ug=fNorm.p1*(Lambda(k).ContourList_px(i,2)-min(Lambda(k).ContourList_px(In,2)))+fNorm.p2;
    %                     ug=fNorm.p1;
                        %ug=fNorm.p1*(Lambda(k).ContourList_px(i,2)-min(Lambda(k).ContourList_px(In,2)))^2+fNorm.p2*(Lambda(k).ContourList_px(i,2)-min(Lambda(k).ContourList_px(In,2)))+fNorm.p3;

                        curvN=[-1 ug];
                        curvN=curvN/max(abs(curvN))/2;
                        if inpolygon(Lambda(k).ContourList_px(i,1)+curvN(1),Lambda(k).ContourList_px(i,2)+curvN(2), Lambda(k).ContourList_px(:,1)', Lambda(k).ContourList_px(:,2)')
                            curvN=[1 -ug];
                            curvN=curvN/max(abs(curvN))/2;
                            if inpolygon(Lambda(k).ContourList_px(i,1)+curvN(1),Lambda(k).ContourList_px(i,2)+curvN(2), Lambda(k).ContourList_px(:,1)', Lambda(k).ContourList_px(:,2)')
                                  curvN=[-1 -1];
                            end
                        end

                    end
                    curvNall=[curvNall; curvN];


                    u=[curvN(1) curvN(2)];
                    v=[WindDir.x WindDir.y];%v=[mean_VXn mean_VYn];
                    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
                    ThetaInDegrees = real(acosd(CosTheta));
                    %[theta,R] = cart2pol(curvN(1),curvN(2));
                    if ThetaInDegrees<=SpreadMax%(theta<=-1.8326 | theta>=0.2618)%theta<=-2.3562 | theta>=0.7854
                        tttt(i)=tttt(i)+1;
                        tttt_Nor(i)=1;
                    end
                else 
                    tttt(i)=tttt(i)+1;
                    tttt_Nor(i)=1;
                end
            end

            for i=1:length(Lambda(k).VX_contour_m_s)
                if length(Lambda(k).VX_contour_m_s)>3
                    %assignin('base','xxx',curvNall);
                    %assignin('base','III',i);
                    if curvNall(i,:)==[-1 -1]
                        aa1=Lambda(k).ContourList_px(:,1)-Lambda(k).ContourList_px(i,1);
                        aa2=Lambda(k).ContourList_px(:,2)-Lambda(k).ContourList_px(i,2);
                        tre=sqrt(aa1.^2+aa2.^2);
                        tre(tre==0)=22;
                        [nNorm,In]=mink(tre,4);
                        curvNall(i,1)=mean(curvNall(In,1));
                        curvNall(i,2)=mean(curvNall(In,2));
                        curvN(1)=mean(curvNall(In,1));
                        curvN(2)=mean(curvNall(In,2));
                        [theta,R] = cart2pol(curvN(1),curvN(2));
                        u=[curvN(1) curvN(2)];
                        v=[WindDir.x WindDir.y];%v=[mean_VXn mean_VYn];
                        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
                        ThetaInDegrees = real(acosd(CosTheta));

                        if ThetaInDegrees<=SpreadMax%(theta<=-1.8326 | theta>=0.2618)
                            tttt(i)=tttt(i)+1;
                            tttt_Nor(i)=1;
                        end
                    end
                end
            end


            Lambda(k).curvN=curvNall;
            Lambda(k).lenEl=lenEl;

            Lambda(k).tttt=tttt;

            % Specific criteria
            Lambda(k).tttt_X_grad=tttt_X_grad;
            Lambda(k).tttt_T_grad=tttt_T_grad;
            Lambda(k).tttt_DIR=tttt_DIR;
            Lambda(k).tttt_Nor=tttt_Nor;
            
            Lambda(k).X_grad=X_grad;
            Lambda(k).T_grad=T_grad;

            Velss=[];
            % Find any velocities not pointing outward of the contour and set
            % them to 0.
            Velss(:,1)=Lambda(k).VX_contour_m_s;
            Velss(:,2)=Lambda(k).VY_contour_m_s;
            Contour_at_dt_x = Lambda(k).ContourList_px(:,1)'+ (Lambda(k).VX_contour_m_s./max(abs(Velss'),[],1)/2);
            Contour_at_dt_y = Lambda(k).ContourList_px(:,2)'+ (Lambda(k).VY_contour_m_s./max(abs(Velss'),[],1)/2);

            [IN]    = inpolygon(Contour_at_dt_x,Contour_at_dt_y, Lambda(k).ContourList_px(:,1)', Lambda(k).ContourList_px(:,2)');
            IN      = double(IN);
            t       = find(IN==1);
            IN(t)   = NaN;
            IN      = IN+1;

            t=find(isnan(IN));
            Lambda(k).tttt(t) = Lambda(k).tttt(t)-5;  
        end
    end
    clear vx_region vy_region pixel_list_local vx_region2 vy_region2 vxC_region2 vyC_region2 vxD_region2 vyD_region2 vx_region3 vy_region3  pixel_list_local_contour IN Contour_at_dt_x Contour_at_dt_y

end

if ~isempty(selected_regions)
    for i=1:length(Lambda)
        mean_VX(i) = nanmean(Lambda(i).VX_contour_m_s);
        mean_VY(i) = nanmean(Lambda(i).VY_contour_m_s);
        Lambda(i).mean_VX=mean_VX(i);
        Lambda(i).mean_VY=mean_VY(i);
        
    end
    
    [theta,R] = cart2pol(mean_VX,mean_VY);
    theta_mean = nanmean(theta);
    
    
    % Note check direction
    %t=find(theta>=theta_mean-pi/2 & theta<=theta_mean+pi/2);
    
    t=1:length(Lambda);
    %Determine if area increases
%     t2=[];
%     bw = im2bw(uint16(B(:,:,1)),th./2^16);
%     
%     [Bb2,Ll2] = bwboundaries(bw,'noholes');
%     % get corresponding region properties
%     stats2 = regionprops(Ll2,'Area','Centroid');
%     lll1=[];
%     lll2=[];
%     for jooj=1:length(stats2)
%         lll1=[lll1 stats2(jooj).Centroid(1)];
%         lll2=[lll2 stats2(jooj).Centroid(2)];
%     end
%     %tttt_Area=zeros(1,length(Lambda(k).VX_contour_m_s));
%     for i=1:length(Lambda)
%         tttt_Area=0;
%         if ~isempty(stats2)
%             
%             
%             [val1,idx] = min(sqrt((Lambda(i).Centroid_px(1)-lll1).^2+(Lambda(i).Centroid_px(2)-lll2).^2));
%             
%             Lambda(i).Areass(1,1)=stats2(idx).Area*dx*dx;
%             Lambda(i).Areass(1,2)=Lambda(i).Area_m2*1.00;
%             if stats2(idx).Area*dx*dx>Lambda(i).Area_m2*1.00
%                 Lambda(i).tttt=Lambda(i).tttt+1;
%                 tttt_Area=tttt_Area+1;
%                 
%             end
%             Lambda(i).tttt_Area=tttt_Area;
%         end
%     end
    
    
    % check for glint
%     filenameA=filenames2(iop).name;
%     filenameB=filenames2(iop+1).name;
%     filename_pre=filenames2(iop-3).name;
%     filename_posle=filenames2(iop+3).name;
%     
%     filename_pre1=filenames2(iop-1).name;
%     filename_posle1=filenames2(iop+1).name;
%     
%     [Apre,Bpre,ALL_x_min1,ALL_x_max1,ALL_y_min1,ALL_y_max1] = get_common_FOV_per_pair_V2(dirin,filenameA,filename_pre,PARAM_dx);
%     [Aposle,Bposle,ALL_x_min2,ALL_x_max2,ALL_y_min2,ALL_y_max2] = get_common_FOV_per_pair_V2(dirin,filenameA,filename_posle,PARAM_dx);
%     
%     [Apre1,Bpre1,ALL_x_min11,ALL_x_max11,ALL_y_min11,ALL_y_max11] = get_common_FOV_per_pair_V2(dirin,filenameA,filename_pre1,PARAM_dx);
%     [Aposle1,Bposle1,ALL_x_min21,ALL_x_max21,ALL_y_min21,ALL_y_max21] = get_common_FOV_per_pair_V2(dirin,filenameA,filename_posle1,PARAM_dx);
%     
%     [aaaa,bbbb,ALL_x_min3,ALL_x_max3,ALL_y_min3,ALL_y_max3] = get_common_FOV_per_pair_V2(dirin,filenameA,filenameB,PARAM_dx);
%     aaaa=A;
%     
%     
%     
%     
%     difX1=round((ALL_x_min3-ALL_x_min1)/PARAM_dx);%/PARAM_dx
%     difX2=round((ALL_x_min3-ALL_x_min2)/PARAM_dx);
%     difY1=-round((ALL_y_max3-ALL_y_max1)/PARAM_dx);
%     difY2=-round((ALL_y_max3-ALL_y_max2)/PARAM_dx);
%     
%     difX11=round((ALL_x_min3-ALL_x_min11)/PARAM_dx);%/PARAM_dx
%     difX21=round((ALL_x_min3-ALL_x_min21)/PARAM_dx);
%     difY11=-round((ALL_y_max3-ALL_y_max11)/PARAM_dx);
%     difY21=-round((ALL_y_max3-ALL_y_max21)/PARAM_dx);
%     for k=1:length(Lambda)
%         suma=0;
%         
%         minX=min(Lambda(k).ContourList_px(:,1));
%         maxX=max(Lambda(k).ContourList_px(:,1));
%         minY=min(Lambda(k).ContourList_px(:,2));
%         maxY=max(Lambda(k).ContourList_px(:,2));
%         
%         if max(Lambda(k).tttt)>3
%             shiftX=mean(Lambda(k).VX_contour_m_s)*dt/dx*3;
%             shiftY=mean(Lambda(k).VY_contour_m_s)*dt/dx*3;
%             
%             shiftX1=mean(Lambda(k).VX_contour_m_s)*dt/dx*1;
%             shiftY1=mean(Lambda(k).VY_contour_m_s)*dt/dx*1;
%             
%             Lambda(k).my=shiftY;
%             Lambda(k).mx=shiftX;
%             if (-shiftX + minX -10 +difX1)>1 & (-shiftY + minY -10 +difY1)>1 & (-shiftX + maxX +10 +difX1)<length(Bpre(1,:)) &  (-shiftY + maxY +10 +difY1)<length(Bpre(:,1))
%                 if (-shiftX1 + minX -10 +difX11)>1 & (-shiftY1 + minY -10 +difY11)>1 & (-shiftX1 + maxX +10 +difX11)<length(Bpre1(1,:)) &  (-shiftY1 + maxY +10 +difY11)<length(Bpre1(:,1))
% 
%                     inddx=floor(minX-shiftX+difX1)-5:ceil(maxX-shiftX+difX1)+5;
%                     inddy=floor(minY-shiftY+difY1)-5:ceil(maxY-shiftY+difY1)+5;
%                     inddx1=floor(minX-shiftX1+difX11)-5:ceil(maxX-shiftX1+difX11)+5;
%                     inddy1=floor(minY-shiftY1+difY11)-5:ceil(maxY-shiftY1+difY11)+5;
%                     if (nansum(nansum(Bpre(inddy,inddx)>Threshold))>nansum(nansum(aaaa(minY:maxY,minX:maxX)>Threshold))*0.6) & (nansum(nansum(Bpre1(inddy1,inddx1)>Threshold))>nansum(nansum(aaaa(minY:maxY,minX:maxX)>Threshold))*0.6)
%                         suma=suma+1;
%                     end
%                 end
%             end
%             
%             if (shiftX + minX -10 +difX2)>1 & (shiftY + minY -10 +difY2)>1 & (shiftX + maxX +10 +difX2)<length(Bposle(1,:)) &  (shiftY + maxY +10 +difY2)<length(Bposle(:,1))
%                 if (shiftX1 + minX -10 +difX21)>1 & (shiftY1 + minY -10 +difY21)>1 & (shiftX1 + maxX +10 +difX21)<length(Bposle1(1,:)) &  (shiftY1 + maxY +10 +difY21)<length(Bposle1(:,1))
% 
%                     inddx=floor(minX+shiftX+difX2)-5:ceil(maxX+shiftX+difX2)+5;
%                     inddy=floor(minY+shiftY+difY2)-5:ceil(maxY+shiftY+difY2)+5;
%                     inddx1=floor(minX+shiftX1+difX21)-5:ceil(maxX+shiftX1+difX21)+5;
%                     inddy1=floor(minY+shiftY1+difY21)-5:ceil(maxY+shiftY1+difY21)+5;
%                     if (nansum(nansum(Bposle(inddy,inddx)>Threshold))>nansum(nansum(aaaa(minY:maxY,minX:maxX)>Threshold))*0.6) & (nansum(nansum(Bposle1(inddy1,inddx1)>Threshold))>nansum(nansum(aaaa(minY:maxY,minX:maxX)>Threshold))*0.6)
%                         suma=suma+1;
%                     end
%                 end
%             end
%             
%         end
%         if suma==0
%             Lambda(k).tttt = Lambda(k).tttt-3;
%         end
%     end    




%     if suma==0
%         Lambda(k).tttt = 0;
%     end
    %t=t2;
    Lambda2=Lambda;
    
else
    
    Lambda2=[];
    Lambda=[];
end
return

      
  
    
    
    
    
    



