function [ Hhat_LS] = functionCellFreeLS( H,nbrOfRealizations,M,K,p,tau_p,Pset)

%LS channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs. Note that the covariance and mean are not given.




%INPUT:
%
%H                    = Matrix with dimension M x nbrOfRealzations x K
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
%Hhat_LS           = Matrix with dimension M x nbrOfRealzations x K
%                     where (m,n,k) is the n^th  realization of LS 
%                     channel estimate between AP m and UE k

%Prepare to store MMSE channel estimates
Hhat_LS = zeros(M,nbrOfRealizations,K);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations) + 1i*randn(M,nbrOfRealizations));

        %Go through all UEs
        for k=1:K
            
            inds=Pset(:,k);
            yp=zeros(M,nbrOfRealizations);
            for z=1:length(inds)
                yp = yp + sqrt(p(inds(z)))*tau_p*H(:,:,inds(z));
            end
            yp = yp + sqrt(tau_p)*Np;
            Hhat_LS(:,:,k)= (1/(sqrt(p(k))*tau_p))*yp;
            
        end

end

