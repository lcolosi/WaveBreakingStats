function Threshold=det_Threshold(N_all,arrayHist,peakPercentage)
% Determines threshold for whitecaps from histogram of images
%[temp,IndeksMin]=max(N_all(5:end));
%IndeksMin=IndeksMin+4;
%IndeksMin=1;
N_CS=cumsum(N_all(200:500));
N_CS=N_CS/nanmax(N_CS);
N_CS=1-N_CS;

%Array=gradient(gradient(log(movmean(N_CS,1))));
%Array=movmean(Array,3);
Array=gradient(movmean(gradient(movmean(log(N_CS),11)),11));

[MaxA,Indeks]=max(Array);
ThresholdInd=find(Array(Indeks:end)<(peakPercentage/100)*MaxA,1);
Threshold=(ThresholdInd+Indeks+200+3)*(arrayHist(2)-arrayHist(1));
% IndeksMax=find(N_CS<0.001,1);
% % The tail can be noisy, remove it
% [MaxA,Indeks]=max(Array(IndeksMin:IndeksMax));
% ThresholdInd=find(Array(Indeks+IndeksMin:end)<(peakPercentage/100)*MaxA,1);
% Threshold=(ThresholdInd+IndeksMin+Indeks-1+1)*(arrayHist(2)-arrayHist(1));

