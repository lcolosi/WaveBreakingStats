function meanIm_2=mean2D(meanIm,winSize)

    %%%%
    % meanIm_2=mean2D(meanIm,winSize)
    %
    % Function for computing a 2D moving average of an array using a square
    % kernel of side length winSize. 
    % 
    %
    %   Parameters
    %   ----------
    %   MeanIm   : NxM array; cannot contain NaNs.   
    %   winSize  : Side length of the square kernel.   
    % 
    %   Returns
    %   -------
    %   meanIm_2 : 2D averaged array of the same size as the input array.
    %  
    % 
    %%%%

    % Initialize kernel
    k = ones(winSize,winSize)/winSize^2;

    % Convolve the mean image with the kernel (the 'same' argument returns
    % the central part of the convolution, which is the same size as tempIm)
    meanIm_2=conv2(meanIm,k,'same'); 

end

%% Developmental code 
% % Initialize a ones array 
% tempIm=ones(size(meanIm));
% 
% % Convolve the ones array with the kernel (same returns the central 
% % part of the convolution, which is the same size as tempIm)
% tempIm=conv2(tempIm,k,'same');
% meanIm_2=conv2(meanIm,k,'same')./tempIm;