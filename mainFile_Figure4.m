%This Matlab script can be used to generate Figure 4 in the article:
%
%Ö. Özdogan, E. Björnson and J. Zhang, "Performance of Cell-Free Massive MIMO with Rician Fading 
%and Phase Shifts," in IEEE Transactions on Wireless Communications.
%doi: 10.1109/TWC.2019.2935434
%
%Download article: https://arxiv.org/abs/1903.07335
%
%This is version 1.0 (Last edited: 2019-11-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

%Empty workspace and close figures
clear
close all
clc

%Define the number of Access Points (APs)
M=100;
%Number of UEs
K=40; 
%Pilot length
tau_p=5;
%Select length of coherence block
tau_c=200;

%The cell size 1km x 1km
cellRange=1000;

%Compute the noise power
noiseFigure = 7;
B=20e6;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
sigma2=db2pow(noiseVariancedBm-30); %(in Watt)
%Uplink transmit power per UE (W)
p=0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
pv=p*ones(1,K);


%Select the number of setups with random AP/UE locations
nbrOfSetups=50;
nbrOfRealizations=1;


%Prepare to save simulation results 
userSE_MMSE= zeros(K,nbrOfSetups,1);
userSE_LMMSE= zeros(K,nbrOfSetups,1);
userSE_LS= zeros(K,nbrOfSetups,1);


%Prepare to save simulation results 
userSE_MMSE_LSFD= zeros(K,nbrOfSetups,1);
userSE_LMMSE_LSFD= zeros(K,nbrOfSetups,1);
userSE_LS_LSFD= zeros(K,nbrOfSetups,1);





for m=1:length(M)
    
    %Deploy APs randomly
    APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
    %For single layer decoding set all large-scale fading coefficients to 1
    A_singleLayer=reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    
    %Go through all setups
    for n=1:nbrOfSetups
       %Deploy UEs and generate the covariance and mean matrices
       [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K,cellRange,APpositions,sigma2,1);
      
       
       %Channel generations for each UE-AP pair
       [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K );
       
       %Pilot allocation 
       [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K,M(m),pv,tau_p);
       
       %Second Layer Decoding
       %Generate Large-scale fading coefficients for phase-aware MMSE,
       %Linear MMSE and LS estimators
       [ A_MMSE,A_LMMSE,A_LS] = functionLSFD( R,HMeanWithoutPhase,M(m),K,pv,tau_p,Pset);
       
       
       %Second layer decoding Theoretical Spectral Efficiency (SE) computation
       %SE with MMSE estimator
       [SE_CC_MMSE_LSFD] = functionTheoreticalCellFreeULSE_MMSE( R,HMeanWithoutPhase,A_MMSE,M(m),K,pv,tau_p,tau_c,Pset);
       %SE with LMMSE estimator
       [SE_CC_LMMSE_LSFD] = functionTheoreticalCellFreeULSE_LMMSE( R,HMeanWithoutPhase,A_LMMSE,M(m),K,pv,tau_p,tau_c,Pset);
       %SE with LS estimator
       [SE_CC_LS_LSFD] = functionTheoreticalCellFreeULSE_LS( R,HMeanWithoutPhase,A_LS,M(m),K,pv,tau_p,tau_c,Pset);
       
      
    
       %Single layer decoding Theoretical Spectral Efficiency (SE) computation
       %SE with MMSE estimator
       [SE_CC_MMSE] = functionTheoreticalCellFreeULSE_MMSE( R,HMeanWithoutPhase,A_singleLayer,M(m),K,pv,tau_p,tau_c,Pset);
       %SE with LMMSE estimator
       [SE_CC_LMMSE] = functionTheoreticalCellFreeULSE_LMMSE( R,HMeanWithoutPhase,A_singleLayer,M(m),K,pv,tau_p,tau_c,Pset);
       %SE with LS estimator
       [SE_CC_LS] = functionTheoreticalCellFreeULSE_LS( R,HMeanWithoutPhase,A_singleLayer,M(m),K,pv,tau_p,tau_c,Pset);

       %Store SEs for all UEs
       userSE_MMSE(:,n) = SE_CC_MMSE(:);
       userSE_LMMSE(:,n) = SE_CC_LMMSE(:);
       userSE_LS(:,n) = SE_CC_LS(:);
       
       
       userSE_MMSE_LSFD(:,n) = SE_CC_MMSE_LSFD(:);
       userSE_LMMSE_LSFD(:,n) = SE_CC_LMMSE_LSFD(:);
       userSE_LS_LSFD(:,n) = SE_CC_LS_LSFD(:);
       
       %Output simulation progress
       disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    end
   
  %Output simulation progress
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end
% Plot the simulation results
CDFnumbers = linspace(0,1,K*nbrOfSetups);
figure(1);
hold on; box on;
plot(sort(userSE_MMSE(:)),CDFnumbers,'r-');
plot(sort(userSE_LMMSE(:)),CDFnumbers,'k-');
plot(sort(userSE_LS(:)),CDFnumbers,'b-.');
plot(sort(userSE_MMSE_LSFD(:)),CDFnumbers,'r--');
plot(sort(userSE_LMMSE_LSFD(:)),CDFnumbers,'k--');
plot(sort(userSE_LS_LSFD(:)),CDFnumbers,'b:');
xlabel('SE per UE [bit/s/Hz]');
ylabel('Cumulative Distribution Function (CDF)');
legend('MMSE estimator','LMMSE estimator','LS estimator','Location','SouthEast');

