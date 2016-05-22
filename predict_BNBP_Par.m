function logprob = predict_BNBP_Par(Xtrain,output,Kdex,IsKnowKall)
% Calculate the predictive likelihood of a count vector under 
% a beta-negative binomial process random count matrix
% Check Equations (33) and (12), Section 3.4 and Appendix D of arXiv:1404.3331
if ~exist('IsKnowKall','var')
    IsKnowKall= false;
end
a0=1e-3;
b0=1e-3;

gamma0=output.gamma0;
c=output.c;
r_i=output.r_i ;
p_k=output.p_k ;
p_star=output.p_star;
n_dot_k=output.n_dot_k;
K=size(Xtrain,1);
nk=zeros(K,1);
nk(Kdex)=n_dot_k;

probsum = sum(-log(1-p_k))+p_star;


r = sum(r_i);
logprob = zeros(1,size(Xtrain,2));

nk_zero_dex = (nk==0);
Kold = length(n_dot_k);
for j=1:size(Xtrain,2)
    Xj = nonzeros(Xtrain(:,j));
    
    % Calculate rj using EM
    lj=max(length(Xj),1);
    rj = (lj+a0-1)/(b0+ probsum);
    for iterEM=1:20
        lj = max(sum(rj*(psi(rj+Xj)-psi(rj))),1);
        rj = (lj+a0-1)/(b0+ probsum);
    end
    
    if IsKnowKall
        %Equation (21) for a finite vocabulary
        temp =  sum(betaNB_pdf_log(full(Xtrain(:,j)),rj,nk+gamma0/K,c+r));
    else
        %Equation (33) and (12) for an infinite vocabulary
        Xj = Xtrain(:,j)>0;
        Xj2 = Xj & nk_zero_dex;
        Knew = nnz(Xj2);
        temp =  sum(betaNB_pdf_log(full(Xtrain(Kdex,j)),rj,n_dot_k,c+r));
        temp = temp + Knew*log(gamma0)-gamma0*(psi(c+r+rj)-psi(c+r));
        if Knew>0
            temp =temp+ gammaln(Kold+1) - gammaln(Kold+Knew+1);
            temp =temp - gammaln(Knew+1);
            
            temp = temp + sum(digamma_pdf_log_unnormalized(Xtrain(Xj2,j),rj,c+r));
            
        %    temp = temp + gammaln(c+r+rj)*Knew- sum(gammaln(c+Xtrain(Xj2,j)+r+rj))...
        %        +sum(gammaln(rj+Xtrain(Xj2,j)))-sum(log(Xtrain(Xj2,j)))-Knew*gammaln(rj);
        end
    end
    logprob(j)=temp;
end
