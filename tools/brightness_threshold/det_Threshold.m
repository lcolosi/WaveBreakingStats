function Threshold = det_Threshold(N_all,arrayHist,peakPercentage)

    %%%%
    % Threshold = det_Threshold(N_all,arrayHist,peakPercentage)
    %
    % Function for determining the threshold for whitecaps from 
    % histogram of images.     
    %
    %   Parameters
    %   ---------- 
    %   N_all          : 
    %   arrayHist      : 
    %   peakPercentage : 
    % 
    %   Returns
    %   ------- 
    %   Threshold      : 
    % 
    %%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute complementary cumulative distribution function 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_CS=cumsum(N_all(200:500));
    N_CS=N_CS/max(N_CS,[],'omitnan');
    N_CS=1-N_CS;
    
    % 
    Array=gradient(movmean(gradient(movmean(log(N_CS),11)),11));
    
    %
    [MaxA,Indeks]=max(Array);
    ThresholdInd=find(Array(Indeks:end)<(peakPercentage/100)*MaxA,1);
    Threshold=(ThresholdInd+Indeks+200+3)*(arrayHist(2)-arrayHist(1));

    % Loop through tracks
    for i=1:length(tracks)
        for j=1:length(N(i).hi(:,1))
            if N(i).hi(j,1)==0
                Threshold(i).th(j)=0;
            else
                N2=(N(i).hi(j,:));
                N_CS=cumsum(N2(200:500));
                N_CS=N_CS/max(N_CS,[],'omitnan');
                N_CS=1-N_CS;
                ind1=find(N_CS<0.5,1,'first');
                [res_x, idx_of_result]=knee_pt(N_CS(ind1:301),ind1:301);
                Threshold(i).th(j)=(199+ind1+idx_of_result)*16;
            end
            
        end
    end
    

end

%% Development Code

%[temp,IndeksMin]=max(N_all(5:end));
%IndeksMin=IndeksMin+4;
%IndeksMin=1;

%Array=gradient(gradient(log(movmean(N_CS,1))));
%Array=movmean(Array,3);

% IndeksMax=find(N_CS<0.001,1);
% % The tail can be noisy, remove it
% [MaxA,Indeks]=max(Array(IndeksMin:IndeksMax));
% ThresholdInd=find(Array(Indeks+IndeksMin:end)<(peakPercentage/100)*MaxA,1);
% Threshold=(ThresholdInd+IndeksMin+Indeks-1+1)*(arrayHist(2)-arrayHist(1));