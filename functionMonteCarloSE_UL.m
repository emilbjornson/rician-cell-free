function [SE_MR] = functionMonteCarloSE_UL(Hhat,H,A,tau_c,tau_p,nbrOfRealizations,M,K,p)

%Computes the uplink SE by Monte Carlo simulations.

%INPUT:
%Hhat                = Matrix with dimension M x nbrOfRealizations x K where
%                      (m,n,k) is the n^th  realization of the channel estimate
%                       between AP m and UE k (when one of the MMSE,LMMSE or LS
%                       estimators is used.)
%H                   = Matrix with dimension M x nbrOfRealzations x K
%                      where (m,n,k) is the n^th channel realization
%                      between AP m and UE k 
%A                   = Diagonal matrix with dimension M x M x K where (:,:,k)
%                      is the LSFD coefficients of UE k (when one of the MMSE
%                      ,LMMSE or LS estimators is used.)
%tau_c               = Length of the coherence block
%tau_p               = Pilot length
%nbrOfRealizations   = Number of realizations
%M                   = Number of APs
%K                   = Number of UEs 
%p                   = 1xK vector, uplink power at each UE
%
%
%OUTPUT:
%
%SE_MR              = Vector with dimension K x 1 where (k) is the SE of UE k
                   



%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store simulation results
SE_MR = zeros(K,1);

Term1x=zeros(K,nbrOfRealizations);
Term2x=zeros(K,nbrOfRealizations);
Term3x=zeros(K,nbrOfRealizations);

    % Go through all channel realizations
    for n = 1:nbrOfRealizations

        %All channels for one realization
        Hallj = reshape(H(:,n,:),[M K]);

        for k = 1:K

            Term1x(k,n) =Hhat(:,n,k)'*A(:,:,k)'*H(:,n,k);
            Term2x(k,n)= sum(p.*abs(Hhat(:,n,k)'*A(:,:,k)'*Hallj).^2) ;
            Term3x(k,n)=norm(A(:,:,k)*Hhat(:,n,k))^2;
        
        end

    end

Term1=real(mean(Term1x,2));
Term2=mean(Term2x,2);
Term3=real(mean(Term3x,2));

%Calculate the SE of each UE k
for k=1:K
    SE_MR(k)= prelogFactor*log2(1+ (p(k)*abs(Term1(k))^2)/(Term2(k) - p(k)*abs(Term1(k))^2 +    Term3(k) )    )  ;
end


end
