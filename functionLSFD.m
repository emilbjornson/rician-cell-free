function [ a_MMSE,a_LMMSE,a_LS] = functionLSFD( R,HMean,M,K,p,tau_p,Pset)

%This function returns the LSFD coefficients in diagonal matrix form
%for two layer decoding.


%INPUT:
%R                   = Matrix with dimension M x M x K
%                      where (:,:,k) is the spatial correlation matrix
%                      between APs and UE k
%                     
%HMean               = Matrix with dimension M x K where (m,k) is the
%                      channel mean between AP m and UE k 
%                      with random phase shifts
%M                    = Number of APs
%K                    = Number of UEs 
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set

%OUTPUT:
%a_MMSE               = Matrix with dimension M x M x K where (:,:,k) 
%                       is the LSFD coefficient of UE k when phase-aware
%                       MMSE estimator is applied.
%                   
%                     
%a_LMMSE               = Matrix with dimension M x M x K where (:,:,k) 
%                       is the LSFD coefficient of UE k when LMMSE estimor
%                       is applied.
%                      
%a_LS                   = Matrix with dimension M x M x K where (:,:,k) 
%                       is the LSFD coefficient of UE k when LMMSE estimor
%                       is applied.




%Prepare to store results
A_MMSE_part1=zeros(M,M,K,K);
A_LMMSE_part1=zeros(M,M,K,K);
A_LS_part1=zeros(M,M,K,K);

%Prepare to store the matrix R' in the paper
Rp=zeros(M,M,K);
Lk=zeros(M,M,K);

eyeM=eye(M);
a_MMSE=zeros(M,M,K);
a_LMMSE=zeros(M,M,K);
a_LS=zeros(M,M,K);

for k=1:K
    Lk(:,:,k)=diag(HMean(:,k).^2);
    Rp(:,:,k)=R(:,:,k)+ Lk(:,:,k);
end

%Go through all UEs
for k = 1:K
    
    %Compute the matrix that is inverted in the MMSE&LMMSE estimators
    inds=Pset(:,k);
    PsiInv=zeros(M,M);
    PsiInv_LMMSE=zeros(M,M);
    %Go through the UEs that use the same pilot
    for z=1:length(inds)
        PsiInv = PsiInv +p(inds(z))*tau_p*R(:,:,inds(z)) ;
        PsiInv_LMMSE = PsiInv_LMMSE +p(inds(z))*tau_p*Rp(:,:,inds(z)) ;
    end
    PsiInv = PsiInv  + eyeM;
    PsiInv_LMMSE = PsiInv_LMMSE  + eyeM;
    
    %Calculate the given terms in the LSFD receiver vector
    %for MMSE&LMMSE
    Omegap=Rp(:,:,k)/PsiInv_LMMSE*Rp(:,:,k);
    Zk=p(k)*tau_p*R(:,:,k)/PsiInv*R(:,:,k)+ Lk(:,:,k);
    bk_LMMSE=p(k)*tau_p*diag(Omegap);
    bk_LS=diag(Rp(:,:,k));
    
    %Go through all UEs
    %Non-coherent interference (i=k')
    for l=1:K
        
        A_MMSE_part1(:,:,k,l)=p(l)*Zk*Rp(:,:,l);
        A_LMMSE_part1(:,:,k,l)=p(l)*p(k)*tau_p*Rp(:,:,l)*Omegap;
        A_LS_part1(:,:,k,l)=p(l)*(1/(p(k)*tau_p))*PsiInv_LMMSE*Rp(:,:,l);
        
        %Coherent interference (If there is pilot contamination)
        %Calculating for the phase-aware MMSE estimator
        if any(l==Pset(:,k)) && l~=k
            %MMSE function LSFD calculation
            A_MMSE_part1(:,:,k,l)= A_MMSE_part1(:,:,k,l)...
                + p(l)*p(l)*p(k)*tau_p*tau_p*diag(R(:,:,l)/PsiInv*R(:,:,k))*diag(R(:,:,l)/PsiInv*R(:,:,k))';
            
        end
        %Check the case when l == k
        if l == k
            A_MMSE_part1(:,:,k,k)=  A_MMSE_part1(:,:,k,k)- p(k)* Lk(:,:,k)*Lk(:,:,k);
            
        end
        
        
        
        %Coherent interference (If there is pilot contamination)
        %Calculating for LMMSE and LS estimators
        if any(l==Pset(:,k))
            
            %LMMSE function LSFD calculation
            A_LMMSE_part1(:,:,k,l)=   A_LMMSE_part1(:,:,k,l) ...
                + p(l)*p(l)*p(k)*tau_p*tau_p*(R(:,:,l)*R(:,:,l)/PsiInv_LMMSE*Omegap...
                +2*Omegap/PsiInv_LMMSE*Lk(:,:,l)*R(:,:,l)...
                +(diag(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))*diag(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))' )...
                -(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))*(Rp(:,:,l)/PsiInv_LMMSE*Rp(:,:,k))'   );
            
            
            %LS function LSFD calculation
            A_LS_part1(:,:,k,l)= A_LS_part1(:,:,k,l) + p(l)*(p(l)/p(k))*(R(:,:,l)*R(:,:,l) + 2*Lk(:,:,l)*R(:,:,l)  ...
                + (diag(Rp(:,:,l))*diag(Rp(:,:,l))') - Rp(:,:,l)*Rp(:,:,l) );
            
            
        end
        
    end
    %Producing LSFD coefficients in vector for phase-aware MMSE, LMMSE and
    %LS estimators
    A_MMSE=sum(A_MMSE_part1(:,:,k,:),4) + Zk;
    A_LMMSE=sum(A_LMMSE_part1(:,:,k,:),4) - p(k)*(bk_LMMSE*bk_LMMSE') +p(k)*tau_p*Omegap;
    A_LS=sum(A_LS_part1(:,:,k,:),4)-p(k)*(bk_LS*bk_LS') + (1/(p(k)*tau_p))*PsiInv_LMMSE;
    
    
    
    %Creating the diagonal matrix form of LSFD coefficients (for each UE k)
    a_MMSE(:,:,k)=diag(A_MMSE\diag(Zk));
    a_LMMSE(:,:,k)=diag(A_LMMSE\bk_LMMSE);
    a_LS(:,:,k)=diag(A_LS\bk_LS);
    
end

end

