function logprob = digamma_log(x,r,c)
logprob = gammaln(x+r)+gammaln(c+r) - log(x) -gammaln(c+r+x)-gammaln(r) -log(psi(c+r)-psi(c));
