function logprob = poisspdf_log(x,lambda)
logprob = x.*log(lambda)-lambda-gammaln(x+1);