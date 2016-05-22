function u=rand_TNB(a,p)
u=1;
PMF=1;
while (1)
    if a~=0
        prob = exp(gammaln(u-a) - sum(log(1:u)))/gamma(-a)*p^u/((1-p)^a-1);    
    else
        prob = -1/log(1-p)*p^u/u;
    end
    if prob/(PMF)>rand(1)        
        break        
    end
    PMF= PMF-prob;
    u=u+1;
  %  [u/10000,prob*1000000/PMF]
end

