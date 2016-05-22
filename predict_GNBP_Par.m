function logprob = predict_GNBP_Par(Xtrain,output,Kdex,LogF,IsKnowKall)
% Calculate the predictive likelihood of a count vector under 
% a gamma-negative binomial process random count matrix
% Check Equations (29) and (12), Section 3.4 and Appendix D of arXiv:1404.3331
if ~exist('IsKnowKall','var')
    IsKnowKall= false;
end

gamma0=output.gamma0;
c=output.c;
p_i=output.p_i ;
r_k=output.r_k ;
r_star=output.r_star;
l_k = output.l_k;

%q = sum(-log(max(1-p_i,realmin)));
q = output.q;

K=size(Xtrain,1);

lk = zeros(K,1);
lk(Kdex) = l_k; %output{i}.n_dot_k;

sumprob = sum(r_k)+r_star;

nk_dex= lk>0;
nk_zero_dex = (lk==0);
Kold = nnz(lk);
parfor j=1:size(Xtrain,2)
    
    % Calculate pj using EM
    pj = 1/(1+ (1*1e-3+sumprob)/(1e-3+sum(Xtrain(:,j))));
    qj = -log(max(1-pj,realmin));
    
    if IsKnowKall
        %Equation (20) for a finite vocabulary
        temp = sum(gammaNB_pdf_log(Xtrain(:,j),lk+gamma0/K,c+q,pj,LogF));
     
    else
        %Equation (29) and (12) for an infinite vocabulary
        Xj = Xtrain(:,j)>0;
        Xj2 = Xj & nk_zero_dex;
        Knew = nnz(Xj2);
        
        temp = sum(gammaNB_pdf_log(Xtrain(Kdex,j),l_k,c+q,pj,LogF));
        temp = temp + Knew*log(gamma0)-gamma0*(log(c+q+qj)-log(c+q));
               %%%%%%
        if Knew>0
            temp = temp + gammaln(Kold+1) - gammaln(Kold+Knew+1);
            temp = temp - gammaln(Knew+1);
            %temp = temp + sum(log(sumStirling(Xtrain(Xj2,j),c+q,pj,LogF)));
            temp = temp + sum(sumStirling_log(Xtrain(Xj2,j),c+q,pj,LogF));
            %temp = temp + gammaln(c+r+rj)*Knew- sum(gammaln(c+Xtrain(Xj2,j)+r+rj))...
            %    +sum(gammaln(rj+Xtrain(Xj2,j)))-sum(log(Xtrain(Xj2,j)))-Knew*gammaln(rj);
        end
    end
    logprob(j)= temp;
end