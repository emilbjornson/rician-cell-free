function [Hhat_LMMSE] = functionCellFreeLMMSE(R,HMeanWithoutPhase,H,nbrOfRealizations,M,K,p,tau_p,Pset)

%LMMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.Note that HMean should be without phase HMeanWithoutPhase




%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMeanWithoutPhase   = M x K channel mean matrix (without random phase shifts)
%                      before introducing random phases
%H                   = Matrix with dimension M x nbrOfRealzations x K
%                     where (m,n,k) is the n^th channel realization
%                     between AP m and UE k (only used for generating y_p)
%nbrOfRealizations    = Number of realizations
%M                    = Number of APs
%K                    = Number of UEs 
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set
%
%
%OUTPUT:
%
%Hhat_LMMSE           = Matrix with dimension M x nbrOfRealzations x K
%                     where (m,n,k) is the n^th  realization of LMMSE 
%                     channelestimate between AP m and UE k


%Prepare to store LMMSE channel estimates
Hhat_LMMSE = zeros(M,nbrOfRealizations,K);

%Store identity matrix of size M x M
eyeM = eye(M);
%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations) + 1i*randn(M,nbrOfRealizations));

%Prepare to store R' matrix
Lk=zeros(M,M,K);
Rp=zeros(M,M,K);
for k=1:K
    Lk(:,:,k)=diag(HMeanWithoutPhase(:,k).^2);
    Rp(:,:,k)=R(:,:,k) + Lk(:,:,k);
end

    for k=1:K
       %Compute the matrix that is inverted in the LMMSE estimator
            inds=Pset(:,k);
            PsiInv_LMMSE=zeros(M,M);
            yp=zeros(M,nbrOfRealizations);
           
            for z=1:length(inds)
                PsiInv_LMMSE = PsiInv_LMMSE +p(inds(z))*tau_p*Rp(:,:,inds(z)) ;
                yp = yp + sqrt(p(inds(z)))*tau_p*H(:,:,inds(z));
               
            end
            %Compute the matrix that is inverted in the LMMSE estimator
            PsiInv_LMMSE = PsiInv_LMMSE  + eyeM;
            yp = yp + sqrt(tau_p)*Np;
        
        
        
        %Save the result
        RPsi = Rp(:,:,k)  / PsiInv_LMMSE;
        Hhat_LMMSE(:,:,k) = sqrt(p(k))*RPsi*yp;


    end
    



