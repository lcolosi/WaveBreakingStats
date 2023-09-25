function meanIm_2=mean2D(meanIm,winSize)

% A utility 2d moving average function.

tempIm=ones(size(meanIm));
tempIm=conv2(tempIm,ones(winSize,winSize)/winSize^2,'same');

meanIm_2=conv2(meanIm,ones(winSize,winSize)/winSize^2,'same')./tempIm;