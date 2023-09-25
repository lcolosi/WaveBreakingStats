% Determine mean brightness and std for corrected images
function [meanOriginal,RM_Nr,meanIm,stdIm] = determineMeanStd(beginDif,endDif,dirRaw,D_Im,arrayHist)

parfor j=beginDif:1:endDif
    
    if mod(j,5000)==0
       j 
    end
    a=imread([dirRaw D_Im(j).name]);
    a=double(im2gray(a));
    a2=imflatfield(a,300);

    
    ImSeries2=sort(a(:));
    ImSeries2=ImSeries2(1:floor(0.8*length(ImSeries2)));
    meanOriginal(j)=mean(ImSeries2);
end

[tempp,RM_Nr]=rmoutliers(meanOriginal(beginDif:endDif));
aTemp=imread([dirRaw D_Im(beginDif).name]);
meanIM2=zeros(size(aTemp));
counter=zeros(size(aTemp));
parfor j=beginDif:1:endDif
    if ~ismember(j,RM_Nr-1+beginDif)
        if mod(j,5000)==0
            j
        end
        a=imread([dirRaw D_Im(j).name]);
        a=double(im2gray(a));
        a2=imflatfield(a,300);
        
        a2(a2>median(a2(:))+5*std(a2(:)) | a2<median(a2(:))-5*std(a2(:)))=nan;
        counter=counter+double(~isnan(a2));
        a(isnan(a2))=0;
        meanIM2=meanIM2+a;
    end
end
meanIm=meanIM2./counter;

stdIM2=zeros(size(aTemp));
counter=zeros(size(aTemp));
parfor j=beginDif:1:endDif
    if ~ismember(j,RM_Nr-1+beginDif)
        if mod(j,5000)==0
            j
        end
        a=imread([dirRaw D_Im(j).name]);
        a=double(im2gray(a));
        a2=imflatfield(a,300);

        a2(a2>median(a2(:))+5*std(a2(:)) | a2<median(a2(:))-5*std(a2(:)))=nan;
        counter=counter+double(~isnan(a2));
        
        addTemp=((a-meanIm).^2);
        addTemp(isnan(a2))=0;
        
        stdIM2=stdIM2+addTemp;
    end
end
stdIm=(stdIM2./counter).^0.5;

