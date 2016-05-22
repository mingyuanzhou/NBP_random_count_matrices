function logprob = nbinpdf_log(x,r,p)
logprob = zeros(size(x));
if length(p)==1
    logprob(x==0) = r(x==0).*log(1-p);
    x1 = x(x~=0);
    r1 = r(x~=0);
    logprob(x~=0) = gammaln(x1+r1)-gammaln(x1+1)-gammaln(r1)+x1*log(p)+r1.*log(1-p);
else
    logprob(x==0) = r(x==0).*log(1-p(x==0));
    x1 = x(x~=0);
    r1 = r(x~=0);
    logprob(x~=0) = gammaln(x1+r1)-gammaln(x1+1)-gammaln(r1)+x1.*log(p(x~=0))+r1.*log(1-p(x~=0));
end