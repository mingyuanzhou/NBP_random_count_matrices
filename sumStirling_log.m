function y = sumStirling_log(x,c,p,LogF)
if nargin<4
    LogF = LogFmatrix(max(x));
end
%% 
% c=10; p=0.9;
% sum(exp(sumStirling_log(1:1000,c,p)))/(log(c-log(1-p))-log(c))
% sum( sumStirling(1:\infty,c,p,LogF)/(ln(c-ln(1-p))-ln(c))) = 1
%pmf = (p^n \sum_{l=0}^n (|s(n,l)|/n!) Gamma(l)/(c-ln(1-p))^l ) / (ln(c-ln(1-p))-ln(c) ) 
y=zeros(size(x));
% for i=1:length(x)
%     n=x(i);
%     y(i) = sum(exp((LogF{n})'+gammaln(1:n)-(1:n).*log(c-log(1-p))+n*log(p) ));
% end

y(x==1)=log(p)-log(c-log(1-p));
y(x==2)= 2*log(p)-log(2) -log(c-log(1-p)) +log(1+ 1./(c-log(1-p)));
y(x==3)= 3*log(p) -log(6) -log(c-log(1-p))  + log(2+ 3./(c-log(1-p)) +2./(c-log(1-p)).^2); 
y(x>3) = sumStirling2_log(x(x>3),c,p,LogF);
end

function y = sumStirling2_log(x,c,p,LogF)
y=zeros(size(x));
for i=1:length(x)
     n=x(i);
     y(i) = log(sum(exp((LogF{n})'+gammaln(1:n)-(1:n).*log(c-log(1-p))+n*log(p) )));
end
end

