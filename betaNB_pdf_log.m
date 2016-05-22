function y = betaNB_pdf_log(x,r,a,b)

%y = gammaln(r+x)-gammaln(r)-gammaln(x+1)+betaln(b+r,a+x)-betaln(a,b);
% tic
% y = gammaln(r+x)-gammaln(r)-gammaln(x+1)...
%     +gammaln(b+r)+gammaln(a+x)-gammaln(a+b+r+x)...
%     +gammaln(a+b)-gammaln(a)-gammaln(b);
%
% sum(y)
% toc
%tic
if length(a)==1
    a = a*ones(size(x));
end
if 1
    y=zeros(size(x));
    x2 = x(x>1);
    a2 = a(x>1);
    
    a1 = a(x==1);
    a0 = a(x==0);
    
    y(x>1) = gammaln(r+x2)-gammaln(r)-gammaln(x2+1)...
        +gammaln(b+r)+gammaln(a2+x2)-gammaln(a2+b+r+x2)...
        +gammaln(a2+b)-gammaln(a2)-gammaln(b);
    
    y(x==1) = log(r)+gammaln(b+r)+log(a1)-gammaln(a1+b+r+1)...
        +gammaln(a1+b)-gammaln(b);
    
    y(x==0) = gammaln(b+r)-gammaln(a0+b+r)...
        +gammaln(a0+b)-gammaln(b);

    %toc
    
    
else
    x2 = x(x>1);
    a2 = a(x>1);
    
    a1 = a(x==1);
    a0 = a(x==0);
    
    y=sum(gammaln(r+x2)-gammaln(r)-gammaln(x2+1)...
        +gammaln(b+r)+gammaln(a2+x2)-gammaln(a2+b+r+x2)...
        +gammaln(a2+b)-gammaln(a2)-gammaln(b));
    
    y=y+sum(log(r)+gammaln(b+r)+log(a1)-gammaln(a1+b+r+1)...
        +gammaln(a1+b)-gammaln(b));
    
    y=y+sum( gammaln(b+r)-gammaln(a0+b+r)...
        +gammaln(a0+b)-gammaln(b));
end