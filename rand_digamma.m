function n=rand_digamma(r,c,maxIter)
n=1;
PMF=1;
count = 0;
while (1)
    count=count+1;
    
    prob = exp(gammaln(n+r)+gammaln(c+r)-gammaln(c+n+r) -gammaln(r))/n/(psi(c+r)-psi(c));
    if prob/(PMF)>rand(1)        
        break        
    end
    PMF= PMF-prob;
    n=n+1;
    if count>maxIter
        break
    end
  %  [u/10000,prob*1000000/PMF]
end

