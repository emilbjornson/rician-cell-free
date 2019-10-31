function [SE_CC] = functionTheoreticalCellFreeDLSE_LMMSE( R,HMean,D,M,K,p,tau_p,tau_c,Pset)

%Theoretical DL SE calculation for LMMSE estimator (non-coherent transmission)

%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMean               = M x K channel mean matrix (without random phase shifts)
%                     
%D                   = Matrix with dimesion M x M x K where (:,:,k) is the
%                      DL power allocation matrix of the UE k
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
%SE_CC              = Vector with dimension K x 1 where (k) is the DL SE of UE k




%Store identity matrix of size M x M
%assuming only uplink transmission
eyeM = eye(M);

prelogFactor = (tau_c-tau_p)/(tau_c);

CCterm1=zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}

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
    PsiInv_LMMSE = p*tau_p*sum(Rp(:,:,Pset(:,k)),3)  + eyeM;
    Omega=(Rp(:,:,k))/PsiInv_LMMSE*(Rp(:,:,k));
    CCterm1(k)=p*tau_p*trace(D(:,:,k)*Omega) ;
   %Non-coherent interference (i=k')
    for i=1:K  
        
        CCterm2_p1(k,i)=trace(D(:,:,i)*Rp(:,:,k));
        PsiInv_LMMSE_l = p*tau_p*sum(Rp(:,:,Pset(:,i)),3)  + eyeM;
        %Coherent interference (If there is pilot contamination)
        if any(i==Pset(:,k))  
            
            
 CCterm2_p1(k,i)=   CCterm2_p1(k,i) + p*tau_p*trace(D(:,:,i)*(R(:,:,k)*R(:,:,k)+2*Lk(:,:,k)*R(:,:,k))/PsiInv_LMMSE_l);
         
            
            
        end
        
        
        
    end
end
     
       CCterm2=sum(CCterm2_p1,2);
       %Calculate the SE of each UE k  
       for k=1:K
           SE_CC(k)= prelogFactor*log2(1+ (CCterm1(k))/(CCterm2(k) - CCterm1(k)+    1 )  )  ;
       end

end

