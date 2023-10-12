function [Threshold,N,Glint]=local_Thresh(peakPercentage,dirProcessed,D_Im,tracks_Im,trackTag,tracks,arrayHist,meanIm,RM_Nr,Glint,Glint_mask,localStep)
% This function determines the brightness theshold (over which a pixel is 
% considered to be whitecap) for each track during a day.
% Input: 
% 

numCores=12;
poolobj = parpool(numCores);
% nizB=[2,10,14];
for i=11:12%5:length(tracks)
%     i=nizB(iiii);
    if trackTag(i).stable==1
        %     imTemp0=mean2D(meanIm(i).im,15);
        %     imTemp=sort(imTemp0(:));
        %     imTemp=imTemp(1:floor(end/2));
        %     imMedian=median(imTemp);
        %     imStd=std(imTemp);
        
        
        
        % Find biggest connected area corresponding to glint using the
        % selected criteria.
        %     BinaryImage=imTemp0>imMedian+imStd*3;
        %     stats = regionprops(BinaryImage, 'Area','PixelIdxList');
        %     [biggest_area,Indeks] = max( [stats.Area] );
        %     Glint(i).im=stats(Indeks).PixelIdxList;
        
        beginDif=tracks_Im(i).Indices(1)+trackTag(i).range(1)-tracks(i).Indices(1);
        endDif=tracks_Im(i).Indices(2)+trackTag(i).range(2)-tracks(i).Indices(2);
        
        dirProcessed2=[dirProcessed 'Track_' num2str(i) '\'];
        filenames=dir([dirProcessed2 '*.tif']);
        
        dirProcessed2_mask=[dirProcessed 'TrackMask_' num2str(i) '\'];
        filenames_mask=dir([dirProcessed2_mask '*.tif']);
        %for j=beginDif:endDif
        initalNr=1+trackTag(i).range(1)-tracks(i).Indices(1)-1;
%         N(i).hi=[];
%         N(i).hi=nan(,length(arrayHist));
        N2=[];
        %% par dole
        parfor j=1+trackTag(i).range(1)-tracks(i).Indices(1):length(filenames)+trackTag(i).range(2)-tracks(i).Indices(2)
            initalNr=1+trackTag(i).range(1)-tracks(i).Indices(1)-1;
            if RM_Nr(i).nr(j-initalNr)==0
                if mod(j,100)==0
                    j
                end
                %a=imread([dirProcessed '\Track_' num2str(i) '\' D_Im(j).name]);
                a=imread([dirProcessed2 filenames(j).name]);
                a=double(im2gray(a));
                
                a_mask=imread([dirProcessed2_mask filenames_mask(j).name]);
                a_mask=double(im2gray(a_mask));
                a_mask=a_mask*nanmax(meanIm(i).im(:))/256/256;
                glintMag=Glint(i).list;
%                 assignin('base','imageTemp',a_mask);
                BinaryImage=a_mask>glintMag;
%                 assignin('base','imageTemp2',BinaryImage);
                statss = regionprops(BinaryImage, 'Area','PixelIdxList');
                [biggest_area,Indeks] = max( [statss.Area] );
                GlintTemp=statss(Indeks).PixelIdxList;
                a(GlintTemp)=nan;
                %assignin('base','imageTemp2',a);
                
                %N(i).hi(j-initalNr,:)=histcounts(a,arrayHist);
                N2(j,:)=histcounts(a,arrayHist);
                tempA=a(:);
                tempA=sort(tempA(~isnan(tempA)));
                tempA=tempA(1:floor(end/2));
                StdJ(j)=std(tempA);
            end
        end
        N(i).hi=movmean(N2,localStep,1,'omitnan');
        
        for j=1+trackTag(i).range(1)-tracks(i).Indices(1):length(filenames)+trackTag(i).range(2)-tracks(i).Indices(2)
            if RM_Nr(i).nr(j-initalNr)==0
                Threshold(i).th(j)=det_Threshold(N(i).hi(j,:),arrayHist,peakPercentage);
%                 Threshold(i).th(j)=5040+2*StdJ(j);
            end
        end
        %assignin('base','test',Threshold);
    end
end
delete(poolobj)
