function AA=removeGlintArea(AA,Glint,a_mask)

A=AA;
glintMag=Glint;
BinaryImage=a_mask>glintMag;
statss = regionprops(BinaryImage, 'Area','PixelIdxList');
[biggest_area,Indeks] = max( [statss.Area] );
GlintTemp=statss(Indeks).PixelIdxList;
AA(GlintTemp)=nan;

% Only take largest connected area
BinaryImage=~isnan(AA);
statss = regionprops(BinaryImage, 'Area','PixelIdxList');
[biggest_area,Indeks] = max( [statss.Area] );
for i=1:length(statss)
   if i~= Indeks
       GlintTemp=statss(i).PixelIdxList;
       AA(GlintTemp)=nan;
   end
end