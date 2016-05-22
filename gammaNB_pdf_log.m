function y = gammaNB_pdf_log(x,a,c,p,LogF)
%Calculate the log likelihood of the gamma-negative binomial distirbution 
%described in Appendix D of arXiv:1404.3331
%a=1;c=1;p=0.5; sum(exp(gammaNB_pdf_log(0:1000,a,c,p)))

if nargin<5
    LogF = LogFmatrix(max(x));
end
% if length(a)==1
%     for i=1:length(x)
%
%         if x(i)==0
%             y(i)=a*log(c)-a*log(c-log(1-p));
%         elseif x(i)==1
%             %y(i)=log(sum(exp( a*log(c)+log(p)+log(a)-(a+1)*log(c-log(1-p)))));
%             y(i)= a*log(c)+log(p)+log(a)-(a+1)*log(c-log(1-p));
%         else
%             y(i)=log(sum(exp( a*log(c)+x(i)*log(p)-gammaln(a) + (LogF{x(i)})'+ gammaln(a+(1:x(i)))-(a+(1:x(i)))*log(c-log(1-p)))));
%         end
%     end
% else
%     for i=1:length(x)
%
%         if x(i)==0
%             y(i)=a(i)*log(c)-a(i)*log(c-log(1-p));
%         elseif x(i)==1
%             %y(i)=log(sum(exp( a(i)*log(c)+log(p)+log(a(i))-(a(i)+1)*log(c-log(1-p)))));
%             y(i)= a(i)*log(c)+log(p)+log(a(i))-(a(i)+1)*log(c-log(1-p));
%         else
%             y(i)=log(sum(exp( a(i)*log(c)+x(i)*log(p)-gammaln(a(i)) + (LogF{x(i)})'+ gammaln(a(i)+(1:x(i)))-(a(i)+(1:x(i)))*log(c-log(1-p)))));
%         end
%     end
% end

if length(a)==1
    a = a*ones(size(x));
end
y=zeros(size(x));
a0 = a(x==0);
a1 = a(x==1);
a2 = a(x==2);
a3 = a(x==3);

y(x==0) = a0*log(c)-a0*log(c-log(1-p));
y(x==1) = a1*log(c)+log(p)+log(a1)-(a1+1)*log(c-log(1-p));
y(x==2) = a2*log(c)+2*log(p)-log(2) + log(a2) -(a2+1)*log(c-log(1-p))  + log( 1 + (a2+1)./(c-log(1-p)));
y(x==3) = a3*log(c)+3*log(p)-log(6) + log(a3) -(a3+1)*log(c-log(1-p))  + log( 2 + 3*(a3+1)./(c-log(1-p)) +  (a3+1).*(a3+2)./(c-log(1-p))^2 );
y(x>3) = gammaNBpdfLog(x(x>3),a(x>3),c,p,LogF);
end

function logpdf= gammaNBpdfLog(x,a,c,p,LogF)
logpdf=zeros(size(x));
for i=1:length(x)
    logpdf(i)=log(sum(exp( a(i)*log(c)+x(i)*log(p)-gammaln(a(i)) + (LogF{x(i)})'+ gammaln(a(i)+(1:x(i)))-(a(i)+(1:x(i)))*log(c-log(1-p)))));
    %logpdf(i)=a(i)*log(c)+x(i)*log(p)-gammaln(a(i)) + log(sum(exp(  (LogF{x(i)})'+ gammaln(a(i)+(1:x(i)))-(a(i)+(1:x(i)))*log(c-log(1-p)))));
end
end
