%This Matlab script can be used to generate Figure 11 in the article:
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
K=[10,40]; 
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
nbrOfSetups=200;
%Select the number of channel realizations per setup
nbrOfRealizations=1;

%Prepare to save simulation results 
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
        %For pilot allocation
        A_singleLayer=reshape(repmat(eye(M(m)),1,K(k)),M(m),M(m),K(k));
        
        %Go through all setups
        for n=1:nbrOfSetups
            
            
            
            %Deploy UEs and generate the covariance and mean matrices
            [R,HMeanWithoutPhase,channelGain] = functionCellFreeSetup( M(m),K(k),cellRange,APpositions,sigma2,1);
            
            % DL power allocation matrix (Heuristic)
            %The power is allocated proportional to the channel quality
            Dx=[];
            D=zeros(M(m),M(m),K(k));
            for s=1:M(m)
                temp=channelGain(s,:)./sum(channelGain(s,:));
                Dx=[Dx;temp];
            end
            for z=1:K(k)
                D(:,:,z)=rho*diag(Dx(:,z));
            end
            pv=p*ones(1,K(k));
            
            
            %Channel generations for each UE-AP pair
            [H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M(m),nbrOfRealizations,K(k) );
            
            %Pilot allocation (same as UL)
            [Pset] = functionPilotAllocation( R,HMeanWithoutPhase,A_singleLayer,K(k),M(m),pv,tau_p);
            
            %Coherent case power allocation matrix
            [ Dk_MMSE,Dk_LMMSE,Dk_LS] = functionDLPowerControl( R,HMean,M(m),K(k),tau_p,p,rho,Pset,Dx );
            
            %COHERENT TRANSMISSION
            %DL SE with MMSE estimator
            [SE_CC_MMSE] = functionTheoreticalCellFreeDLSE_MMSE_coherent( R,HMeanWithoutPhase,Dk_MMSE,M(m),K(k),p,tau_p,tau_c,Pset);
            %DL SE with LMMSE estimator
            [SE_CC_LMMSE] = functionTheoreticalCellFreeDLSE_LMMSE_coherent( R,HMeanWithoutPhase,Dk_LMMSE,M(m),K(k),p,tau_p,tau_c,Pset);
            %DL SE with LS estimator
            [SE_CC_LS] = functionTheoreticalCellFreeDLSE_LS_coherent( R,HMeanWithoutPhase,Dk_LS,M(m),K(k),p,tau_p,tau_c,Pset);
            
            
            
            if k==1
                %Store SEs for all UEs
                MMSE1(:,n) = SE_CC_MMSE(:);
                LMMSE1(:,n) = SE_CC_LMMSE(:);
                LS1(:,n) = SE_CC_LS(:);
            end
            if k==2
                %Store SEs for all UEs
                MMSE2(:,n) = SE_CC_MMSE(:);
                LMMSE2(:,n) = SE_CC_LMMSE(:);
                LS2(:,n) = SE_CC_LS(:);
                
            end
            
            %Output simulation progress
            disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        end
        
    end
    
    %Output simulation progress
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end


%Plot the simulation results
CDFnumbers = linspace(0,1,K(1)*nbrOfSetups);
CDFnumbers2 = linspace(0,1,K(2)*nbrOfSetups);

figure(1);
hold on; box on;
plot(sort(MMSE1(:)),CDFnumbers,'r-');
plot(sort(LMMSE1(:)),CDFnumbers,'k-');
plot(sort(LS1(:)),CDFnumbers,'b-.');
plot(sort(MMSE2(:)),CDFnumbers2,'r--');
plot(sort(LMMSE2(:)),CDFnumbers2,'k--');
plot(sort(LS2(:)),CDFnumbers2,'b:');
c1=plot(0,0,'b-');
c2=plot(0,0,'b--');
c3=plot(0,0,'b:');
xlabel('SE per UE [bit/s/Hz]');
ylabel('Cumulative Distribution Function (CDF)');
legend('MMSE estimator','LMMSE estimator','LS estimator','Location','SouthEast');
legend([c1 c2 c3],'MMSE','LMMSE','LS','Location','SouthEast');    
 