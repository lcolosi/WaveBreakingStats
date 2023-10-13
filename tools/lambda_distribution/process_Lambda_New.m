function process_Lambda_New(i1,i2,filenames2,dirin,diroutL,Threshold,a,dxx,dtt,WindDir,SpreadMax,Glint,dirinIM_mask,meanIm,RM_Nr,GlobalOrLocal,trackTag,tracks)
%Written by Teodor Vrecica a code that calls the function that determines
%lambda distribution in pair of succesive images, and saves the results.
% Input - 
numCores=14;
poolobj = parpool(numCores);
addpath('utilities')
% parfor i=i1:i2
parfor ijoj=i1:i1+floor((i2-i1)/12)  
     i=i1+(ijoj-i1)*12;
   % i=ijoj;
    if RM_Nr(i-(i1-4))~=1
        % Specify the location of images, and the output folder
        
        PARAM_dt = dtt;
        PARAM_dx = dxx;
        
        % Names of the images for which Lambda is calculated
        filenameA=filenames2(i).name;
        filenameB=filenames2(i+1).name;
        filenameC=filenames2(i+2).name;
        filenameD=filenames2(i+3).name;
        
        [A,B,C,D,ALL_x_min,ALL_x_max,ALL_y_min,ALL_y_max] = get_common_FOV_per_pair_V3(dirin,filenameA,filenameB,filenameC,filenameD,PARAM_dx);

        [a_mask,b_mask,c_mask,d_mask,ALL_x_min,ALL_x_max,ALL_y_min,ALL_y_max] = get_common_FOV_per_pair_V3(dirinIM_mask,filenameA,filenameB,filenameC,filenameD,PARAM_dx);
        a_mask=double(im2gray(a_mask));
        a_mask=a_mask*nanmax(meanIm(:))/256/256;
        % Determine area - overlapping areas of images, with glint removed and
        % with padding of 5m around the edges.
        %Area=sum(sum(A>258))*dxx^2;
        
        A=double(A);
        A=removeGlintArea(A,Glint,a_mask);
        %assignin('base','AA',A);
        
        B=double(B);
        C=double(C);
        D=double(D);
        [Area,AA,BB,CC,DD]=det_area(A,B,C,D,dxx,Glint,a_mask);
        A(isnan(AA))=0;
        B(isnan(BB))=0;
        C(isnan(CC))=0;
        D(isnan(DD))=0;
        
        %% increase brightness
%         testIm=ones(size(A));
%         testIm1=linspace(0.96,1.04,length(A(1,:)));
%         testIm3=linspace(1.02,0.98,length(A(:,1)));
%         
%         
%         A=A.*testIm1.*testIm3;
%         B=B.*testIm1.*testIm3;
%         C=C.*testIm1.*testIm3;
%         D=D.*testIm1.*testIm3;
        
%         A=rot90(A,2);
%         B=rot90(B,2);
%         C=rot90(C,2);
%         D=rot90(D,2);
        % Build gradient images and add them as layers to the original
        % image. The gradient images are set to above 0 (to accomodate the
        % optical flow algorithm.
%          assignin('base','AAA',A);
%          assignin('base','BB',B)
%         assignin('base','CC',C)
%         assignin('base','DD',D)
        
        [fx1,fy1]=gradient(A);
        [fx2,fy2]=gradient(B);
        [fx3,fy3]=gradient(C);
        [fx4,fy4]=gradient(D);
        minYY1=min(fy1(:));
        minYY2=min(fy2(:));
        minYY3=min(fy3(:));
        minYY4=min(fy4(:));
        minYY=min([minYY1 minYY2 minYY3 minYY4]);
        minXX1=min(fx1(:));
        minXX2=min(fx2(:));
        minXX3=min(fx3(:));
        minXX4=min(fx4(:));
        minXX=min([minXX1 minXX2 minXX3 minXX4]);
        A3=cat(3,A,(fx1-1*minXX),(fy1-1*minYY));
        B3=cat(3,B,(fx2-1*minXX),(fy2-1*minYY));
        C3=cat(3,C,(fx3-1*minXX),(fy3-1*minYY));
        D3=cat(3,D,(fx4-1*minXX),(fy4-1*minYY));

        
        
        %%%
        
        %     SpreadMax=90;
        %     WindDir=90;
        %%%
        if GlobalOrLocal==1
            Threshold2=Threshold;
            Thresh_t=(Threshold-4500)*0.2;
            Thresh_x=(Threshold-4500)*0.2;
        else
            Threshold2=Threshold(i);
            Thresh_t=(Threshold(i)-4500)*0.2;
            Thresh_x=(Threshold(i)-4500)*0.2;
        end
        tic
%         [Lambda_Tot,Lambda,Area_Tot] = process_pair_images_New_V2(A3,B3,C3,D3,PARAM_dt,PARAM_dx,a,Threshold2,Thresh_t,Thresh_x,WindDir,SpreadMax,i,dirin,filenames2);
        [Lambda_Tot,Lambda,Area_Tot] = process_pair_images_Inbuilt_V2(A3,B3,C3,D3,PARAM_dt,PARAM_dx,a,Threshold2,Thresh_t,Thresh_x,WindDir,SpreadMax,i,dirin,filenames2);

        toc
        nule=num2str(zeros(1,4-length(num2str(i))));
        nule=strrep(nule, ' ', '');
        saveLambdaDir([diroutL,'LAMBDA_' nule num2str(i),'.mat'],Lambda_Tot,Area);
    end
end

delete(poolobj)

%    
return
