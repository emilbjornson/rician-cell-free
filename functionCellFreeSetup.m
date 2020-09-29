function [ R,HMean,channelGain] = functionCellFreeSetup( M,K,cellRange,APpositions, sigma2, correlatedShadowing)

%The file that generates the cell free setup for given number of APs and
%UEs where the locations of APs and cell size are given.


%INPUT:
%M                    = Number of APs
%K                    = Number of UEs 
%cellRange            = The size of network
%APposiitons          = The locations of APs
%sigma2               = The noise power (in Watt)
%correlatedShadowing  = If set to 1 correlated shadow fading is generated, 
%                        else uncorrelated shadow fading
%
%OUTPUT:
%
%R                   = Matrix with dimension M x M x K
%                     where (:,:,k) is the spatial correlation matrix
%                     between APs and UE k, normalized by noise power
%                     
%HMean               = Matrix with dimension M x K where (m,k) is the
%                      channel mean between AP m and UE k, normalized by
%                      noise power and without random phase shifts
%channelGain         = Matrix with dimension M x K where (m,k) is the total
%                      channel gain between AP m and UE k, normalized by
%                      noise power 


%The minimum allowed distance (Access Point to User Equipment)
dmin=10; 
%The maximum distance of LoS occurence
%maxDistanceLoS=300; 
%The height of Access Point
APheigth=12.5;  
%The height of User Equipment
UEheigth=1.5;  

%The height difference between AP and UE 
verticalDistance=APheigth-UEheigth;


%Dropping all UEs while minimum distance requirement is satisfied.
droppedUEs=0;
%Preparing to save the distances and UE positions
distanceAll=zeros(M,K);
UEpositions=zeros(K,1);


%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-cellRange 0 cellRange],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[M 1]);
%Dropping all UEs
while droppedUEs <K
    
    UEposition=cellRange*(rand(1,1) + 1i*rand(1,1));
    horizontalDistance=sqrt(real(APpositions-UEposition).^2+imag(APpositions-UEposition).^2);
    distance=sqrt(verticalDistance.^2 + horizontalDistance.^2);
    
    if isempty(distance(distance<dmin))
        droppedUEs=droppedUEs+1;
        distanceAll(:,droppedUEs)=distance;
        
        %Store UE positions
        UEpositions(droppedUEs)=UEposition;
    end
    
end


%Calculating the distances between all AP and UE pairs
for k = 1:K
    distanceAll(:,k) = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
end

%Prepare to create correlated shadow fading 
%Prepare to save AP-AP and UE-UE distances
distanceAP=zeros(M,M); 
distanceUE=zeros(K,K);  
%Distances between APs
 for m=1:M
    distanceAP(:,m)=sqrt( real(APpositions - APpositions(m)).^2 + imag(APpositions - APpositions(m)).^2         ); 
 end
%Distances between UEs
 for k=1:K
    distanceUE(:,k)=sqrt( real(UEpositions - UEpositions(k)).^2 + imag(UEpositions - UEpositions(k)).^2         ); 
 end

 %Decorrelation distance
 decorr=100; 
 %Creating covarince functions for each pair
 covMatrixAP=2.^(-distanceAP/decorr);
 covMatrixUE=2.^(-distanceUE/decorr);
 %Creating shadow fading realizations
 shadowFadingAP=sqrtm(covMatrixAP)*randn(M,1); 
 shadowFadingUE=sqrtm(covMatrixUE)*randn(K,1);
 %Shadow fading parameter
 delta=0.5;
 %Create the resulting shadow fading matrix
 shadowFadingMatrix=zeros(M,K);
 for k=1:K
     for m=1:M
        shadowFadingMatrix(m,k)=sqrt(delta)*shadowFadingAP(m) +sqrt(1-delta)*shadowFadingUE(k);
     end
 end
 
 
 %Scale with variance 
 sigma_sf=8;
 if correlatedShadowing==1
    %Correlated Shadow Fading Matrix
    shadowFading=10.^(sigma_sf*shadowFadingMatrix/10); 
 else
     %Uncorrelated Shadow Fading Matrix
     shadowFading=10.^(sigma_sf*randn(M,K)/10); 
 end
 
 
 
 
%Path-loss calculation "Cost 231 Walfish-Ikegami Model"
%before adding shadow fading
betaLoS=30.18 +26*log10(distanceAll);
betaNLoS=34.53+ 38*log10(distanceAll);


%probLoS=(rand<((maxDistanceLoS-distanceAll)./maxDistanceLoS));
%Each pair has a LoS
probLoS=ones(size(distanceAll));
%Calculate the distance based Rician Factor 
ricianFactor=10.^(1.3-0.003*distanceAll);

%Prepare to save channel gain in dB
channelGaindB=zeros(M,K);
channelGain_LoS=zeros(M,K);
channelGain_NLoS=zeros(M,K);

%Save the channel gains (in this setup each pair has a LoS path)
channelGaindB(probLoS==1)=-betaLoS(probLoS==1);
channelGaindB(probLoS==0)=-betaNLoS(probLoS==0);
channelGain=(db2pow(channelGaindB).*shadowFading)/sigma2;

%Scale with Rician factor
channelGain_LoS(probLoS==1)= sqrt(ricianFactor(probLoS==1)./(ricianFactor(probLoS==1) +1 )).*sqrt(channelGain(probLoS==1));
channelGain_NLoS(probLoS==1)=(1./(ricianFactor(probLoS==1) +1 )).*(channelGain(probLoS==1));
channelGain_NLoS(probLoS==0)=channelGain(probLoS==0); %note that probLoS is always one in the manuscript
            
%Prepare to save  channel covariance and LoS matrix
R=zeros(M,M,K);
HMean=zeros(M,K);
%HMeanWithoutPhase=zeros(M,K);

for k=1:K
    R(:,:,k)=diag(channelGain_NLoS(:,k));
    HMean(:,k)=channelGain_LoS(:,k); %note that sqrt is already taken in the prev. step.
end


end

