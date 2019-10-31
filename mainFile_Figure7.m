%This Matlab script can be used to generate Figure 7 in the article:
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


%Define the range of number of Access Points (APs)
M=40:20:100;
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
%Downlink transmit power (initial value)
rho=p;

%Select the number of setups with random AP/UE locations
nbrOfSetups=25;
%Select the number of channel realizations per setup
nbrOfRealizations=100;

%Prepare to save simulation results for Theoretical and Monte-Carlo
%MR for Monte-Carlo and CC for theoretical
%NC for Non-Coherent
sumSE_MR_MMSE= zeros(length(M),nbrOfSetups);
sumSE_CC_MMSE= zeros(length(M),nbrOfSetups);

sumSE_MR_LS= zeros(length(M),nbrOfSetups);
sumSE_CC_LS= zeros(length(M),nbrOfSetups);
sumSE_MR_LMMSE= zeros(length(M),nbrOfSetups);
sumSE_CC_LMMSE= zeros(length(M),nbrOfSetups);

sumSE_MR_MMSE_NC= zeros(length(M),nbrOfSetups);
sumSE_CC_MMSE_NC= zeros(length(M),nbrOfSetups);

sumSE_MR_LS_NC= zeros(length(M),nbrOfSetups);
sumSE_CC_LS_NC= zeros(length(M),nbrOfSetups);
sumSE_MR_LMMSE_NC= zeros(length(M),nbrOfSetups);
sumSE_CC_LMMSE_NC= zeros(length(M),nbrOfSetups);



%Go through the number of APs
for m=1:length(M)
    
    %Deploy APs randomly
    APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
    %For pilot allocation
    A_singleLayer=reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    
    %Go through all setups
    for n=1:nbrOfSetups
        
        %Prepare to store the power allocation matrix
        D=zeros(M(m),M(m),K);
        %Deploy UEs and generate the covariance and mean matrices
        [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K,cellRange,APpositions,sigma2,1);
        
        %Non-coherent case: DL power allocation matrix (Heuristic)
        %The power is allocated proportional to the channel quality
        Dx=[];
        for s=1:M(m)
            temp=channelGain(s,:)./sum(channelGain(s,:));
            Dx=[Dx;temp];
        end
        for k=1:K
            D(:,:,k)=rho*diag(Dx(:,k));
        end
        pv=p*ones(1,K);
        
        %Create channel generations for each UE-AP pair
        [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K );
        %Pilot allocation (same as UL)
        [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K,M(m),pv,tau_p);
        
        %Channel estimation 
        %Channel estimation with phase-aware MMSE estimator
        [Hhat_MMSE] =functionCellFreeMMSE(R,HMean,H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
        %Channel estimation with LMMSE estimator
        [Hhat_LMMSE] =functionCellFreeLMMSE(R,HMeanWithoutPhase,H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
        %Channel estimation with LS estimator
        [Hhat_LS] = functionCellFreeLS( H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
        
        %COHERENT TRANSMISSION
        %Coherent case power allocation matrix
        [ Dk_MMSE,Dk_LMMSE,Dk_LS] = functionDLPowerControl( R,HMean,M(m),K,tau_p,p,rho,Pset ,Dx);
        
        %Coherent Transmission
        %DL SE with MMSE estimator
        [SE_MR_MMSE] = functionMonteCarloSE_DL_coherent(Hhat_MMSE,Dk_MMSE,H,tau_c,tau_p,nbrOfRealizations,M(m),K);
        [SE_CC_MMSE] = functionTheoreticalCellFreeDLSE_MMSE_coherent( R,HMeanWithoutPhase,Dk_MMSE,M(m),K,p,tau_p,tau_c,Pset);
        %DL SE with LMMSE estimator 
        [SE_MR_LMMSE] = functionMonteCarloSE_DL_coherent(Hhat_LMMSE,Dk_LMMSE,H,tau_c,tau_p,nbrOfRealizations,M(m),K);
        [SE_CC_LMMSE] = functionTheoreticalCellFreeDLSE_LMMSE_coherent( R,HMeanWithoutPhase,Dk_LMMSE,M(m),K,p,tau_p,tau_c,Pset);
        %UL SE with LS estimator 
        [SE_MR_LS] = functionMonteCarloSE_DL_coherent(Hhat_LS,Dk_LS,H,tau_c,tau_p,nbrOfRealizations,M(m),K);
        [SE_CC_LS] = functionTheoreticalCellFreeDLSE_LS_coherent( R,HMeanWithoutPhase,Dk_LS,M(m),K,p,tau_p,tau_c,Pset);
      
        %Average of all users
        sumSE_CC_MMSE(m,n) = mean(SE_CC_MMSE,1);
        sumSE_MR_MMSE(m,n) = mean(SE_MR_MMSE,1);
        
        sumSE_CC_LMMSE(m,n) = mean(SE_CC_LMMSE,1);
        sumSE_MR_LMMSE(m,n) = mean(SE_MR_LMMSE,1);
        
        sumSE_CC_LS(m,n) = mean(SE_CC_LS,1);
        sumSE_MR_LS(m,n) = mean(SE_MR_LS,1);

        
        %NON-COHERENT TRANSMISSION
        %DL SE with MMSE estimator 
        [SE_MR_MMSE_NC] = functionMonteCarloSE_DL_noncoherent(Hhat_MMSE,H,D,tau_c,tau_p,M(m),K);
        [SE_CC_MMSE_NC] = functionTheoreticalCellFreeDLSE_MMSE( R,HMeanWithoutPhase,D,M(m),K,p,tau_p,tau_c,Pset);
        
        %DL SE with LMMSE estimator 
        [SE_MR_LMMSE_NC] = functionMonteCarloSE_DL_noncoherent(Hhat_LMMSE,H,D,tau_c,tau_p,M(m),K);
        [SE_CC_LMMSE_NC] = functionTheoreticalCellFreeDLSE_LMMSE( R,HMeanWithoutPhase,D,M(m),K,p,tau_p,tau_c,Pset);
        %DL SE with LS estimator 
        [SE_MR_LS_NC] = functionMonteCarloSE_DL_noncoherent(Hhat_LS,H,D,tau_c,tau_p,M(m),K);
       
      
        %Average of all users
        sumSE_CC_MMSE_NC(m,n) = mean(SE_CC_MMSE_NC,1);
        sumSE_MR_MMSE_NC(m,n) = mean(SE_MR_MMSE_NC,1);
        
        sumSE_CC_LMMSE_NC(m,n) = mean(SE_CC_LMMSE_NC,1);
        sumSE_MR_LMMSE_NC(m,n) = mean(SE_MR_LMMSE_NC,1);
        
      
        sumSE_MR_LS_NC(m,n) = mean(SE_MR_LS_NC,1);
        
        
        
        
        
        
        
        %Output simulation progress
         disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        
        
    end
    %Output simulation progress
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
        
end


%Plot the simulation results
figure;
c1=plot(M,mean(sumSE_MR_MMSE,2),'ks','MarkerSize',5);
hold on
c2=plot(M,mean(sumSE_CC_MMSE,2),'r->');
hold on
plot(M,mean(sumSE_MR_LMMSE,2),'ks','MarkerSize',5);
hold on
c3=plot(M,mean(sumSE_CC_LMMSE,2),'r-x');
hold on
plot(M,mean(sumSE_MR_LS,2),'ks','MarkerSize',4)
hold on
c4=plot(M,mean(sumSE_CC_LS,2),'r-o');
hold on
plot(M,mean(sumSE_MR_MMSE_NC,2),'ks','MarkerSize',5);
hold on
c5=plot(M,mean(sumSE_CC_MMSE_NC,2),'k-->');
hold on
plot(M,mean(sumSE_MR_LMMSE_NC,2),'ks','MarkerSize',5);
hold on
c6=plot(M,mean(sumSE_CC_LMMSE_NC,2),'k--x');
hold on
plot(M,mean(sumSE_MR_LS_NC,2),'ks','MarkerSize',4)
hold on
plot(M,mean(sumSE_CC_LMMSE_NC,2),'k--o');
 hold on
plot(40,0,'k--');
 hold on
% c4=plot(40,0,'k-');
xlabel('Number of APs, M');
ylabel('Average DL SE [bit/s/Hz]');
 legend([c2 c3  c4 c5 c6 c1 ],'Coherent trans. with MMSE','Coherent trans. with LMMSE','Coherent trans. with LS',...
     'Non-coherent trans. with MMSE','Non-coherent trans. with LMMSE/LS', 'Monte-Carlo','Location','southeast');
    
    
    
    
    