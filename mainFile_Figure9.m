%This Matlab script can be used to generate Figure 9 in the article:
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
tau_p=[5,20];
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
pv=p*ones(1,K);
%Downlink transmit power (initial value)
rho=p;

%Select the number of setups with random AP/UE locations
nbrOfSetups=200;
%Select the number of channel realizations per setup
nbrOfRealizations=1;


%Prepare to save simulation results 
userSE_MMSE= zeros(K,nbrOfSetups,length(tau_p));
userSE_LMMSE= zeros(K,nbrOfSetups,length(tau_p));
userSE_LS= zeros(K,nbrOfSetups,length(tau_p));
%Prepare to save simulation results 
userSE_MMSE_NC= zeros(K,nbrOfSetups,length(tau_p));
userSE_LMMSELS_NC= zeros(K,nbrOfSetups,length(tau_p));









for m=1:length(M)
    
    %Deploy APs randomly
    APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
    %For pilot allocation
    A_singleLayer=reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    
    %Go through all setups
    for n=1:nbrOfSetups
        %Go through all pilot lengths
        for t=1:length(tau_p)
            
            
            
            %Deploy UEs and generate the covariance and mean matrices
            [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K,cellRange,APpositions,sigma2,1);
            
            %Non-coherent case: DL power allocation matrix (Heuristic)
            %The power is allocated proportional to the channel quality
            D=zeros(M(m),M(m),K);
            Dx=[];
            
            for s=1:M(m)
                temp=channelGain(s,:)./sum(channelGain(s,:));
                Dx=[Dx;temp];
            end
            for k=1:K
                D(:,:,k)=rho*diag(Dx(:,k));
            end
            
            %Create channel generations for each UE-AP pair
            [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K );
            %Pilot allocation (same as UL)
            [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K,M(m),pv,tau_p(t));
            
            %COHERENT TRANSMISSION
            %Coherent case power allocation matrix
            [ Dk_MMSE,Dk_LMMSE,Dk_LS] = functionDLPowerControl( R,HMean,M(m),K,tau_p(t),p,rho,Pset,Dx );
            
            %Coherent Transmission
            %DL SE with MMSE estimator
            [SE_CC_MMSE] = functionTheoreticalCellFreeDLSE_MMSE_coherent( R,HMeanWithoutPhase,Dk_MMSE,M(m),K,p,tau_p(t),tau_c,Pset);
            %DL SE with LMMSE estimator
            [SE_CC_LMMSE] = functionTheoreticalCellFreeDLSE_LMMSE_coherent( R,HMeanWithoutPhase,Dk_LMMSE,M(m),K,p,tau_p(t),tau_c,Pset);
            %DL SE with LS estimator
            [SE_CC_LS] = functionTheoreticalCellFreeDLSE_LS_coherent( R,HMeanWithoutPhase,Dk_LS,M(m),K,p,tau_p(t),tau_c,Pset);
            
            
            
            %Store SEs for all UEs
            userSE_MMSE(:,n,t) = SE_CC_MMSE(:);
            userSE_LMMSE(:,n,t) = SE_CC_LMMSE(:);
            userSE_LS(:,n,t) = SE_CC_LS(:);
            
            
            
        end
        
        
        %Output simulation progress
        disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        
        
    end
    %Output simulation progress
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
    
end



MMSE1= userSE_MMSE(:,:,1);
LMMSE1= userSE_LMMSE(:,:,1);
LS1= userSE_LS(:,:,1) ;


MMSE2= userSE_MMSE(:,:,2);
LMMSE2= userSE_LMMSE(:,:,2);
LS2= userSE_LS(:,:,2) ;

% Plot the simulation results
CDFnumbers = linspace(0,1,K*nbrOfSetups);
figure(1);
hold on; box on;
plot(sort(MMSE1(:)),CDFnumbers,'r-');
plot(sort(LMMSE1(:)),CDFnumbers,'r--');
plot(sort(LS1(:)),CDFnumbers,'r:');
plot(sort(MMSE2(:)),CDFnumbers,'k-');
plot(sort(LMMSE2(:)),CDFnumbers,'k--');
plot(sort(LS2(:)),CDFnumbers,'k:');
c1=plot(0,0,'b-');
c2=plot(0,0,'b--');
c3=plot(0,0,'b:');
xlabel('SE per UE [bit/s/Hz]');
ylabel('Cumulative Distribution Function (CDF)');
legend([c1 c2 c3],'MMSE','LMMSE','LS','Location','SouthEast');    
    
    
    
    