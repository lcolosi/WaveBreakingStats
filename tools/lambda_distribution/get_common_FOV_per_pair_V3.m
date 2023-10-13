function [A,B,C,D,ALL_x_min,ALL_x_max,ALL_y_min,ALL_y_max] = get_common_FOV_per_pair_V3(dirin,filename1,filename2,filename3,filename4,PARAM)

[a1,r1]=geotiffread([dirin filename1]);
%a1=rgb2gray(a1);
xmin1=(r1.XWorldLimits(1));
xmax1=(r1.XWorldLimits(2));
ymin1=(r1.YWorldLimits(1));
ymax1=(r1.YWorldLimits(2));

[a2,r2]=geotiffread([dirin filename2]);
%a2=rgb2gray(a2);
xmin2=(r2.XWorldLimits(1));
xmax2=(r2.XWorldLimits(2));
ymin2=(r2.YWorldLimits(1));
ymax2=(r2.YWorldLimits(2));

[a3,r3]=geotiffread([dirin filename3]);
%a2=rgb2gray(a2);
xmin3=(r3.XWorldLimits(1));
xmax3=(r3.XWorldLimits(2));
ymin3=(r3.YWorldLimits(1));
ymax3=(r3.YWorldLimits(2));

[a4,r4]=geotiffread([dirin filename4]);
%a2=rgb2gray(a2);
xmin4=(r4.XWorldLimits(1));
xmax4=(r4.XWorldLimits(2));
ymin4=(r4.YWorldLimits(1));
ymax4=(r4.YWorldLimits(2));

% overlapping area

ALL_x_min =max([xmin1, xmin2,xmin3, xmin4]);
ALL_x_max =min([xmax1, xmax2,xmax3, xmax4]);
ALL_y_min =max([ymin1, ymin2,ymin3, ymin4]);
ALL_y_max =min([ymax1, ymax2,ymax3, ymax4]);

% find corresponding area for each image
x1_b=round((ALL_x_min-xmin1)/PARAM);
x1_e=round((-ALL_x_max+xmax1)/PARAM);
x2_b=round((ALL_x_min-xmin2)/PARAM);
x2_e=round((-ALL_x_max+xmax2)/PARAM);

x3_b=round((ALL_x_min-xmin3)/PARAM);
x3_e=round((-ALL_x_max+xmax3)/PARAM);
x4_b=round((ALL_x_min-xmin4)/PARAM);
x4_e=round((-ALL_x_max+xmax4)/PARAM);

y1_e=round((ALL_y_min-ymin1)/PARAM);
y1_b=round((-ALL_y_max+ymax1)/PARAM);
y2_e=round((ALL_y_min-ymin2)/PARAM);
y2_b=round((-ALL_y_max+ymax2)/PARAM);

y3_e=round((ALL_y_min-ymin3)/PARAM);
y3_b=round((-ALL_y_max+ymax3)/PARAM);
y4_e=round((ALL_y_min-ymin4)/PARAM);
y4_b=round((-ALL_y_max+ymax4)/PARAM);


A = a1(1+y1_b:end-y1_e,1+x1_b:end-x1_e);
B = a2(1+y2_b:end-y2_e,1+x2_b:end-x2_e);
C = a3(1+y3_b:end-y3_e,1+x3_b:end-x3_e);
D = a4(1+y4_b:end-y4_e,1+x4_b:end-x4_e);


