function [ Dk_MMSE,Dk_LMMSE,Dk_LS] = functionDLPowerControl( R,HMean,M,K,tau_p,p,rho,Pset,Dx )
    %Coherent case: DL power control coefficients based on MMSE, LMMSE and LS estimates
    
    
    %Prepare to store results
    Zk_vect=zeros(M,K);
    Omegap_vect=zeros(M,K);
    eyeM=eye(M);
    Lk=zeros(M,M,K);
    Rp=zeros(M,M,K);
    Rp_vect=zeros(M,K);
    Dk_MMSE=zeros(M,M,K);
    Dk_LMMSE=zeros(M,M,K);
    Dk_LS=zeros(M,M,K);
    %Preparation
    for k=1:K
        Lk(:,:,k)=diag(abs(HMean(:,k)).^2);
        Rp(:,:,k)=R(:,:,k) + Lk(:,:,k);
        Rp_vect(:,k)=diag(Rp(:,:,k));
    end
    
    
    %Go through all UEs
    for k = 1:K
        PsiInv = (p*tau_p*sum(R(:,:,Pset(:,k)),3) + eyeM);
        PsiInv_LMMSE = p*tau_p*sum(Rp(:,:,Pset(:,k)),3)  + eyeM;
        Zk_vect(:,k)=diag(p*tau_p*R(:,:,k)/PsiInv*R(:,:,k)+ Lk(:,:,k));
        Omegap_vect(:,k)=p*tau_p*diag((Rp(:,:,k))/PsiInv_LMMSE*(Rp(:,:,k)));
        
    end
    %Calculate the power allocation matrices for MMSE, LMMSE and LS
    %estimators
    for k=1:K
    Dk_MMSE(:,:,k)=diag(Dx(:,k)*rho./Zk_vect(:,k));
    Dk_LMMSE(:,:,k)=diag(Dx(:,k)*rho./Omegap_vect(:,k));
    Dk_LS(:,:,k)=diag(Dx(:,k)*rho./Rp_vect(:,k));
    
    end
    
end

