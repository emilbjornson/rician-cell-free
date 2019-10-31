function [SE_CC] = functionTheoreticalCellFreeULSE_MMSE( R,HMean,A,M,K,p,tau_p,tau_c,Pset)
        
%Theoretical SE calculation for phase-aware MMSE estimator

%INPUT:
%R                   = M x M x K channel covariance matrix 
%HMean               = M x K channel mean matrix (with random phase shifts)
%                     
%A                   = Diagonal matrix with dimension M x M x K where (:,:,k)
%                      is the LSFD coefficients of UE k (when  MMSE
%                      estimator is used.)
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
prelogFactor=(tau_c -tau_p)/tau_c;


CCterm1=zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm3=zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
%CCterm2=zeros(K,1); %Store E{|v^H_k h_k'|^2}
CCterm2_p1=zeros(K,K);
SE_CC=zeros(K,1);%Store the result
Lk=zeros(M,M,K);
Omega=zeros(M,M,K);
for k=1:K
     Lk(:,:,k)=diag(abs(HMean(:,k)).^2);
end


    %Go through all UEs
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE estimator
        inds=Pset(:,k);
        PsiInv=zeros(M,M);
        for z=1:length(inds)
            PsiInv = PsiInv +p(inds(z))*tau_p*R(:,:,inds(z)) ;
        end
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = PsiInv  + eyeM;
        
        
        Omega(:,:,k)=R(:,:,k)/PsiInv*R(:,:,k);
        Zk=p(k)*tau_p*R(:,:,k)/PsiInv*R(:,:,k)+ Lk(:,:,k);
        CCterm1(k)=trace(A(:,:,k)*Zk);
        CCterm3(k)=trace(A(:,:,k)'*A(:,:,k)*Zk);
        %Non-coherent interference (i=k')
        for l=1:K  
           

                CCterm2_p1(k,l)=p(l)*trace(A(:,:,k)'*Zk*(R(:,:,l)+Lk(:,:,l))*A(:,:,k));
                %Coherent interference (If there is pilot contamination)
                if any(l==Pset(:,k)) && l~=k 

                    CCterm2_p1(k,l)= CCterm2_p1(k,l)+ p(l)*p(l)*p(k)*tau_p*tau_p*abs(trace(A(:,:,k)*R(:,:,l)/PsiInv*R(:,:,k)))^2;
                
                end
               if l == k
                    CCterm2_p1(k,l)= CCterm2_p1(k,l) -  p(k)*trace(A(:,:,k)'*Lk(:,:,k)*Lk(:,:,k)*A(:,:,k));
                
                
                end
        

        end



    end

CCterm2=sum(CCterm2_p1,2);
%Calculate the SE of each UE k    
for k=1:K
    SE_CC(k)= prelogFactor*log2(1+ (p(k)*abs(CCterm1(k))^2)/(abs(CCterm2(k)) +  CCterm3(k) )  )  ;
end


end

