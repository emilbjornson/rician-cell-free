function [ SE_CC] = functionTheoreticalCellFreeDLSE_LS_coherent( R,HMeanWithoutPhase,Dk,M,K,p,tau_p,tau_c,Pset)

%Theoretical DL SE calculation for LS estimator (Coherent transmission)

%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMean               = M x K channel mean matrix (without random phase shifts)
%                     
%Dk                  = Matrix with dimesion M x M x K where (:,:,k) is the
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

CCterm1=zeros(K,1);

%Store E{v^H_k h_k}= E{|v_k|^2}
%CCterm2=zeros(K,1); %Store E{|v^H_k h_k'|^2}
CCterm2_p1=zeros(K,K);
SE_CC=zeros(K,1);%Store the result
Lk=zeros(M,M,K);
Rp=zeros(M,M,K);

for k=1:K
    Lk(:,:,k)=diag(HMeanWithoutPhase(:,k).^2);
    Rp(:,:,k)=R(:,:,k)+ Lk(:,:,k);
end

%Go through all UEs
for k = 1:K
    Dkm=sqrtm(Dk(:,:,k));
    %Compute the matrix that is inverted in the LMMSE estimator
    %PsiInv_LS= p*tau_p*sum(Rp(:,:,Pset(:,k)),3)  + eyeM;
    CCterm1(k)=trace(Dkm*Rp(:,:,k)) ;
    
    %Non-coherent interference (l=k')
    for l=1:K  
        PsiInv_LS_l= p*tau_p*sum(Rp(:,:,Pset(:,l)),3)  + eyeM;
        CCterm2_p1(k,l)=(1/(p*tau_p))*trace(Dk(:,:,l)*PsiInv_LS_l*Rp(:,:,k));
        
        
        if any(l==Pset(:,k))  %Coherent interference (If there is pilot contamination)
            
            CCterm2_p1(k,l)= CCterm2_p1(k,l)+trace(Dk(:,:,l)*(R(:,:,k)*R(:,:,k) + 2*Lk(:,:,k)*R(:,:,k)))  ...
                + trace(sqrtm(Dk(:,:,l))*Rp(:,:,k))^2 - trace(Dk(:,:,l)*Rp(:,:,k)*Rp(:,:,k)) ;
           
        end
        
    end
    
    
    
end

     
       CCterm2=sum(CCterm2_p1,2);
      %Calculate the SE of each UE k   
       for k=1:K
           SE_CC(k)= prelogFactor*log2(1+ (p*abs(CCterm1(k))^2)/(p*CCterm2(k) - p*abs(CCterm1(k))^2 +    1 )    )  ;
       end

end

