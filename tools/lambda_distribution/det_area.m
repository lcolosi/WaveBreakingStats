function [Area,AA,BB,CC,DD]=det_area(A,B,C,D,dxx,Glint,a_mask)
% Determine area used for estimating breaking distributions. First
% overlapping part of the images is determined, the area of high glint is
% cropped out, and boundaries (5m around the edges) are cropped out.

A(A==0)=nan;
B(B==0)=nan;
C(C==0)=nan;
D(D==0)=nan;
A(isnan(B))=nan;
A(isnan(C))=nan;
A(isnan(D))=nan;


glintMag=Glint;
%assignin('base','imageTemp',a_mask);
BinaryImage=a_mask>glintMag;
statss = regionprops(BinaryImage, 'Area','PixelIdxList');
[biggest_area,Indeks] = max( [statss.Area] );
GlintTemp=statss(Indeks).PixelIdxList;
A(GlintTemp)=nan;

% Only take largest connected area
BinaryImage=~isnan(A);
statss = regionprops(BinaryImage, 'Area','PixelIdxList');
[biggest_area,Indeks] = max( [statss.Area] );
for i=1:length(statss)
   if i~= Indeks
       GlintTemp=statss(i).PixelIdxList;
       A(GlintTemp)=nan;
   end
end

B(isnan(A))=nan;
C(isnan(A))=nan;
D(isnan(A))=nan;
AA=A;
BB=B;
CC=C;
DD=D;


lengthW=floor(5/dxx)*2;
A=conv2(A,ones(lengthW,lengthW),'same');
%assignin('base','AA',~isnan(A));
Area=sum(~isnan(A),'all')*dxx^2;
