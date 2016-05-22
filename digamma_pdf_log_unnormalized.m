function y = digamma_pdf_log_unnormalized(x,r,c)

%exp(y)/psi(r+c)-psi(c) is the pmf for digamma
%r=1;c=1; sum(exp(digamma_pdf_log_unnormalized(1:100000,r,c)))/(psi(r+c)-psi(c))
y=zeros(size(x));

y(x==1) = log(r)-log(c+r);
y(x==2) = log(r)+log(r+1)-log(c+r)-log(c+r+1)-log(2);
y(x==3) = log(r)+log(r+1)+log(r+2)-log(c+r)-log(c+r+1)-log(c+r+2)-log(3);
y(x>3) = gammaln(x(x>3)+r)+gammaln(c+r)-gammaln(r)-gammaln(x(x>3)+c+r) - log(x(x>3));    

%