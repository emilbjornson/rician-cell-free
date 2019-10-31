function [coherentx,nonCoherentx] = functionMMSE_interferenceLevels( R,HMean,A,M,K,p,tau_p,Pset)
        
%Check the levels of coherent and non-coherent interference levels (used
%for pilot allocation)


%Store identity matrix of size M x M
eyeM = eye(M);
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

        for l=1:K  %Non-coherent interference (i=k')
           

                nonCoherent(k,l)=p(l)*trace(A(:,:,k)'*Zk*(R(:,:,l)+Lk(:,:,l))*A(:,:,k));

                if any(l==Pset(:,k))  %Coherent interference (If there is pilot contamination)

                    coherent(k,l)=  p(l)*p(l)*p(k)*tau_p*tau_p*abs(trace(A(:,:,k)*R(:,:,l)/PsiInv*R(:,:,k)))^2;
                
                end
 

        end



    end

       coherentx=sum(coherent,2);
       nonCoherentx=sum(nonCoherent,2);
      

end

