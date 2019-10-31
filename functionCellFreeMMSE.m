function [Hhat_MMSE] = functionCellFreeMMSE(R,HMean,H,nbrOfRealizations,M,K,p,tau_p,Pset)

%MMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.


%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMean               = M x K channel mean matrix (with random phase shifts)
%                     
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
%OUTPUT:
%
%Hhat_MMSE            = Matrix with dimension M x nbrOfRealzations x K
%                     where (m,n,k) is the n^th  realization of phase-aware
%                     MMSe channel estimate between AP m and UE k
%                     




%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfRealizations,K);

%Store identity matrix of size M x M
eyeM = eye(M);
%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations) + 1i*randn(M,nbrOfRealizations));

    
  
        for k=1:K
            
            %Compute the matrix that is inverted in the MMSE estimator
            inds=Pset(:,k);
            PsiInv=zeros(M,M);
            yp=zeros(M,nbrOfRealizations);
            yMean=zeros(M,nbrOfRealizations);
            for z=1:length(inds)
                PsiInv = PsiInv +p(inds(z))*tau_p*R(:,:,inds(z)) ;
                yp = yp + sqrt(p(inds(z)))*tau_p*H(:,:,inds(z));
                yMean = yMean + sqrt(p(inds(z)))*tau_p*HMean(:,:,inds(z));
            end
            %Compute the matrix that is inverted in the MMSE estimator
            PsiInv = PsiInv  + eyeM;
            yp = yp + sqrt(tau_p)*Np;
          
            %Save the result
            RPsi = R(:,:,k) / PsiInv;
            Hhat_MMSE(:,:,k) = HMean(:,:,k) + sqrt(p(k))*RPsi*(yp-yMean);
            
        end
       
    end



