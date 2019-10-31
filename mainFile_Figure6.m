%This Matlab script can be used to generate Figure 6 in the article:
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
K=[10,40]; % Only can take {tau_p, 2*tau_p} for this setup.
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


%Select the number of setups with random AP/UE locations
nbrOfSetups=200;
nbrOfRealizations=1;


%Prepare to store simulation results 
MMSE1 = zeros(K(1),nbrOfSetups);
LMMSE1 = zeros(K(1),nbrOfSetups);
LS1 = zeros(K(1),nbrOfSetups);
MMSE2 = zeros(K(2),nbrOfSetups);
LMMSE2 = zeros(K(2),nbrOfSetups);
LS2 = zeros(K(2),nbrOfSetups);





for m=1:length(M)
    
    %Go through all UEs
    for k=1:length(K)
        %Deploy APs randomly
        APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
        %For single layer decoding set all large-scale fading coefficients to 1
        A_singleLayer=reshape(repmat(eye(M(m)),1,K(k)),M(m),M(m),K(k));
        %Create the power vector for all UEs (The uplink power is the same
        %(p)at each UE)
        pv=p*ones(1,K(k));
        
        %Go through all setups
        for n=1:nbrOfSetups
            for t=1:length(tau_p)
                %Deploy UEs and generate the covariance and mean matrices
                [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K(k),cellRange,APpositions,sigma2,1);
               
                %Channel generations for each UE-AP pair
                [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K(k) );
                
                %Pilot allocation 
                [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K(k),M(m),pv,tau_p);
                
                %Second Layer Decoding
                %Generate Large-scale fading coefficients for phase-aware MMSE,
                %Linear MMSE and LS estimators
                [ A_MMSE,A_LMMSE,A_LS] = functionLSFD( R,HMeanWithoutPhase,M(m),K(k),pv,tau_p,Pset);
                
                
                %Second layer decoding Theoretical Spectral Efficiency (SE) computation
                %SE with MMSE estimator
                [SE_CC_MMSE_LSFD] = functionTheoreticalCellFreeULSE_MMSE( R,HMeanWithoutPhase,A_MMSE,M(m),K(k),pv,tau_p,tau_c,Pset);
                %SE with LMMSE estimator
                [SE_CC_LMMSE_LSFD] = functionTheoreticalCellFreeULSE_LMMSE( R,HMeanWithoutPhase,A_LMMSE,M(m),K(k),pv,tau_p,tau_c,Pset);
                %%SE with LS estimator
                %[SE_CC_LS_LSFD] = functionTheoreticalCellFreeULSE_LS( R,HMeanWithoutPhase,A_LS,M(m),K(k),pv,tau_p,tau_c,Pset);
                
                
                
                if k==1
                    %Store SEs for all UEs
                    MMSE1(:,n) = SE_CC_MMSE_LSFD(:);
                    LMMSE1(:,n) = SE_CC_LMMSE_LSFD(:);
                    %                 LS1(:,n) = SE_CC_LS_LSFD(:);
                end
                if k==2
                    %Store SEs for all UEs
                    MMSE2(:,n) = SE_CC_MMSE_LSFD(:);
                    LMMSE2(:,n) = SE_CC_LMMSE_LSFD(:);
                    %                 LS2(:,n) = SE_CC_LS_LSFD(:);
                    
                end
                
                
            end
            %Output simulation progress
            disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        end
    end
    %Output simulation progress
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end


% Plot the simulation results
CDFnumbers = linspace(0,1,K(1)*nbrOfSetups);
CDFnumbers2 = linspace(0,1,K(2)*nbrOfSetups);
figure(1);
hold on; box on;
plot(sort(MMSE1(:)),CDFnumbers,'r-');
plot(sort(LMMSE1(:)),CDFnumbers,'r--');
plot(sort(MMSE2(:)),CDFnumbers2,'k-');
plot(sort(LMMSE2(:)),CDFnumbers2,'k--');
c1=plot(0,0,'b-');
c2=plot(0,0,'b--');
xlabel('SE per UE [bit/s/Hz]');
ylabel('Cumulative Distribution Function (CDF)');
legend([c1 c2],'MMSE  ','LMMSE/LS ','Location','SouthEast');
