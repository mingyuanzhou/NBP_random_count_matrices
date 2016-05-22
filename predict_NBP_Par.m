function logprob = predict_NBP_Par(Xtrain,n_dot_k,gamma0,J,c,IsKnowKall)
% Calculate the predictive likelihood of a count vector under 
% a negative binomial process random count matrix
% Check Equation (7) and (12) of arXiv:1404.3331
if ~exist('IsKnowKall','var')
    IsKnowKall= false;
end
n_dot_k_dex = n_dot_k>0;
n_dot_k_zero_dex = ~n_dot_k_dex;
n_dot_k_nnz = n_dot_k(n_dot_k_dex);
pj=1/(J+c+1);
K=size(Xtrain,1);
Kold = length(n_dot_k_nnz );
logprob=zeros(1,size(Xtrain,2));
parfor j=1:size(Xtrain,2)
    Xj2 = Xtrain(:,j)>0 &(n_dot_k_zero_dex);
    Knew = nnz(Xj2);
    if IsKnowKall
        %Equation (19) for a finite vocabulary
        temp = sum(nbinpdf_log(Xtrain(:,j),n_dot_k+gamma0/K,pj));
    else
        %Equation (7) and (12) for an infinite vocabulary
        temp  = poisspdf_log(Knew,-gamma0*log(max(1-pj,realmin)));
        if Kold>0
           temp  = temp + sum(nbinpdf_log(Xtrain(n_dot_k_dex,j),n_dot_k_nnz,pj));
        end
        if Knew>0
            temp  = temp  + gammaln(Kold+1) - gammaln(Knew+Kold+1);
            temp  = temp + sum(Logrithmic_pdf_log(Xtrain(Xj2,j),pj));
        end
    end
    logprob(j) = temp;
    %logprob(j) = predict_NBP(Xtrain(:,j),n_kCi,gamma0Ci,JCi,cCi,n_kCi_NNZ, n_kCi_NNZ_dex);
end

