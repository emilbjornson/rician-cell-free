function [SE_MR] = functionMonteCarloSE_DL_noncoherent(Hhat,H,D,tau_c,tau_p,M,K)

%Computes the DL SE by Monte Carlo simulations for coherent transmission.

%INPUT:
%Hhat                = Matrix with dimension M x nbrOfRealizations x K where
%                      (m,n,k) is the n^th  realization of the channel estimate
%                       between AP m and UE k (when one of the MMSE,LMMSE or LS
%                       estimators is used.)
%H                   = Matrix with dimension M x nbrOfRealzations x K
%                      where (m,n,k) is the n^th channel realization
%                      between AP m and UE k 
%D                   = Matrix with dimesion M x M x K where (:,:,k) is the
%                      DL power allocation matrix of the UE k
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
rho=0.1;


%Prepare to store simulation results for interference powers

SE_MR=zeros(K,1);
Term1xx=zeros(M,K);
Term1=zeros(K,1);


hhat_mean=zeros(M,K);

for m=1:M
    for k=1:K
    hhat_mean(m,k)=mean(abs(Hhat(m,:,k)).^2);
    end
    
end


%All estimates are normalized 
for k=1:K
    for m=1:M
    Hhat(m,:,k)=Hhat(m,:,k)./sqrt(hhat_mean(m,k));
    end
end



    

for k = 1:K

    for m=1:M
       
        Term1xx(m,k)=D(m,m,k)*abs(mean(conj(H(m,:,k)).*Hhat(m,:,k))).^2;
        for i=1:K
        Term2xx(m,k,i)=D(m,m,i)*(mean(abs(conj(H(m,:,k)).*Hhat(m,:,i)).^2));
        end
        
    end
    Term1(k)=rho*sum(Term1xx(:,k));
   
end

 Term2=rho*sum(sum(Term2xx,3),1);


%Calculate the SE of each UE k

    for k=1:K

        SE_MR(k)= prelogFactor*log2( 1 + Term1(k)/(Term2(k) -Term1(k)+ 1)  ) ;

    end



end




