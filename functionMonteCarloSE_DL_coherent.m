function [SE_MR] = functionMonteCarloSE_DL_coherent(Hhat,Dk,H,tau_c,tau_p,nbrOfRealizations,M,K)

%Computes the DL SE by Monte Carlo simulations for coherent transmission.

%INPUT:
%Hhat                = Matrix with dimension M x nbrOfRealizations x K where
%                      (m,n,k) is the n^th  realization of the channel estimate
%                       between AP m and UE k (when one of the MMSE,LMMSE or LS
%                       estimators is used.)
%Dk                  = Matrix with dimesion M x M x K where (:,:,k) is the
%                      DL power allocation matrix of the UE k
%H                   = Matrix with dimension M x nbrOfRealzations x K
%                      where (m,n,k) is the n^th channel realization
%                      between AP m and UE k 
%tau_c               = Length of the coherence block
%tau_p               = Pilot length
%nbrOfRealizations   = Number of realizations
%M                   = Number of APs
%K                   = Number of UEs 

%
%
%OUTPUT:
%
%SE_MR              = Vector with dimension K x 1 where (k) is the DL SE of UE k










%Compute the prelog factor assuming only downlink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);



%Prepare to store simulation results for interference powers
SE_MR=zeros(K,1);
Term1x=zeros(K,nbrOfRealizations);
Term2x=zeros(K,K,nbrOfRealizations);
Term2mean=zeros(K,K);


%Go through all realizations
for n = 1:nbrOfRealizations
    
    %All channels for one realization
    Hhatall = reshape(Hhat(:,n,:),[M K]);
    
    %Go through all UEs
    for k = 1:K
        Dkm=sqrtm(Dk(:,:,k));
        Term1x(k,n) =Hhat(:,n,k)'*Dkm*H(:,n,k);
        
        for i=1:K
            Term2x(k,i,n)= abs(Hhatall(:,i)'*sqrtm(Dk(:,:,i))*H(:,n,k)).^2 ;
        end
        
   
    end
    
end

Term1=abs(mean(Term1x,2)).^2;


for k=1:K
    for i=1:K
   Term2mean(k,i)=mean(Term2x(k,i,:));
    end
end

Term2=sum(Term2mean,2);
%Calculate the SE of each UE k
    for k=1:K

        SE_MR(k)= prelogFactor*log2( 1 + Term1(k)/(Term2(k)- Term1(k)+ 1)  ) ;

    end



end




