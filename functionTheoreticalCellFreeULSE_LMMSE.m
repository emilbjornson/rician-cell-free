function [SE_CC] = functionTheoreticalCellFreeULSE_LMMSE( R,HMean,A,M,K,p,tau_p,tau_c,Pset)

%Theoretical UL SE calculation for LMMSE estimator

%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMean               = M x K channel mean matrix (without random phase shifts)
%                     
%A                   = Diagonal matrix with dimension M x M x K where (:,:,k)
%                      is the LSFD coefficients of UE k (when LMMSE estimator is used.)
%M                   = Number of APs
%K                   = Number of UEs 
%p                   = 1xK vector, uplink power at each UE
%tau_p               = Pilot length
%tau_c               = Length of the coherence block
%Pset                = Pilot allocation set
%
%
%OUTPUT:
%
%SE_CC              = Vector with dimension K x 1 where (k) is the SE of UE k

%Store identity matrix of size M x M
%assuming only uplink transmission
eyeM = eye(M);

prelogFactor = (tau_c-tau_p)/(tau_c);

CCterm1=zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm3=zeros(K,1);
%CCterm2=zeros(K,1); %Store E{|v^H_k h_k'|^2}
CCterm2_p1=zeros(K,K);
SE_CC=zeros(K,1);%Store the result
Lk=zeros(M,M,K);
Rp=zeros(M,M,K);

for k=1:K
    Lk(:,:,k)=diag(HMean(:,k).^2);
    Rp(:,:,k)=R(:,:,k) + Lk(:,:,k);
end

%Go through all UEs
for k = 1:K
    
    
    %Compute the matrix that is inverted in the LMMSE estimator
    inds=Pset(:,k);
    PsiInv_LMMSE=zeros(M,M);
    for z=1:length(inds)
        PsiInv_LMMSE = PsiInv_LMMSE + p(inds(z))*tau_p*Rp(:,:,inds(z)) ;
    end
    %Compute the matrix that is inverted in the LMMSE estimator
    PsiInv_LMMSE = PsiInv_LMMSE  + eyeM;
    
    
    Omega=(Rp(:,:,k))/PsiInv_LMMSE*(Rp(:,:,k));
    
    CCterm1(k)=p(k)*tau_p*trace(A(:,:,k)*Omega) ;
    CCterm3(k)=p(k)*tau_p*trace(A(:,:,k)'*Omega*A(:,:,k));
    %Non-coherent interference (l=k')
    for l=1:K
        
        CCterm2_p1(k,l)=p(l)*p(k)*tau_p*trace(A(:,:,k)'*Rp(:,:,l)*Omega*A(:,:,k));
        
        %Coherent interference (If there is pilot contamination)
        if any(l==Pset(:,k))
         
            CCterm2_p1(k,l) =   CCterm2_p1(k,l) ...
                + p(l)*p(l)*p(k)*tau_p*tau_p*(trace(A(:,:,k)'*R(:,:,l)*R(:,:,l)/PsiInv_LMMSE*Omega*A(:,:,k))...
                +2*trace(A(:,:,k)'*Omega/PsiInv_LMMSE*Lk(:,:,l)*R(:,:,l)*A(:,:,k))...
                +trace(A(:,:,k)*Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))^2 ...
                -trace(A(:,:,k)'*(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))*(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))*A(:,:,k)));
            
            
            
        end
        
        
        
    end
end

CCterm2=sum(CCterm2_p1,2);
%Calculate the SE of each UE k    
for k=1:K
    SE_CC(k)= prelogFactor*log2(1+ (p(k)*abs(CCterm1(k))^2)/(CCterm2(k) - p(k)*abs(CCterm1(k))^2 +    CCterm3(k) )  )  ;
end

end

