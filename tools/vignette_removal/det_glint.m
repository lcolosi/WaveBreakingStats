function [Glint,Glint_mask]=det_glint(meanIm,stdMag,tracks,trackTag)

for i=1:length(tracks)
    i
    if trackTag(i).stable==1
        imTemp=sort(meanIm(i).im(:));
        imTemp=imTemp(1:floor(end/1));
        imTemp = imTemp(~isnan(imTemp));
        imMedian=median(imTemp);
        imStd=std(imTemp);
        
        % Find biggest connected area corresponding to glint using the
        % selected criteria.
%         BinaryImage=meanIm(i).im>imMedian+imStd*stdMag;
%         statss = regionprops(BinaryImage, 'Area','PixelIdxList');
%         [biggest_area,Indeks] = max( [statss.Area] );
%         Glint(i).list=statss(Indeks).PixelIdxList;
        Glint(i).list=imStd*stdMag+imMedian;
        % Same, but for surrounding area

        Glint_mask(i).list=imStd*(stdMag+1)+imMedian;
    end
    
end