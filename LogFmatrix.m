function LogF = LogFmatrix(count_max)
%Recursively calculate S(m,l), Stirling Numbers of the first kind in log scale
%g(m,l)=log(S(m,l))- log(m!)
%Appendix E of arXiv:1404.3331
LogF = cell(count_max,1);
LogF{1}= 0;
for m=2:count_max
    LogF{m} = zeros(m,1);
    LogF{m}(1) = log((m-1)/m) + LogF{m-1}(1);
    LogF{m}(2:m-1) = log((m-1)/m) +  LogF{m-1}(2:m-1) + log(1+exp(LogF{m-1}(1:m-2) - LogF{m-1}(2:m-1) - log(m-1)));
    LogF{m}(m) = LogF{m-1}(m-1) - log(m);    
end