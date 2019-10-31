function [ H,HMean] = functionChannelGeneration( R,HMeanWithoutPhase,M,nbrOfRealizations,K)

%The file that generates the channel realizations for given covariance and
%mean matrices. The outputs are channel realizations and channel means with
%random phase shifts at a coherence block.


%INPUT:
%R                    = M x M x K channel covariance matrix 
%HMeanWithoutPhase    = M x K channel mean matrix 
%                       before introducing random phases
%M                    = Number of APs
%nbrOfRealizations    = Number of realizations
%K                    = Number of UEs 
%
%OUTPUT:
%
%H                   = Matrix with dimension M x nbrOfRealzations x K
%                     where (m,n,k) is the n^th channel realization
%                     between AP m and UE k
%                     
%HMean               = Matrix with dimension M x nbrOfRealzations x K where
%                      (m,n,k) is the n^th realization of the channel mean
%                      between AP m and UE k

%Prepare to store each channel realization   
H=zeros(M,nbrOfRealizations,K); 
W = (randn(M,nbrOfRealizations,K)+1i*randn(M,nbrOfRealizations,K));


%Same HMean and R for all realizations (at each setup) but phases are different at each
%coherence block
HMean=zeros(M,nbrOfRealizations,K); 
HMeanx=reshape(repmat(HMeanWithoutPhase,nbrOfRealizations,1),M,nbrOfRealizations,K); 

%Create phase shifts and store them in a matrix
%uniformly distributed random variables
angles= -pi + 2*pi*rand(M,nbrOfRealizations,K);
phaseMatrix=exp(1i*angles);
  
%Go through all UEs and apply the channel gains to the spatial
%correlation and mean matrices and introduce the phase shifts 
    for k = 1:K
        
        HMean(:,:,k)= phaseMatrix(:,:,k).*HMeanx(:,:,k);
        Rsqrt = sqrtm(R(:,:,k));
        H(:,:,k) = sqrt(0.5)*Rsqrt*W(:,:,k) + HMean(:,:,k);
       
    end
    
end

