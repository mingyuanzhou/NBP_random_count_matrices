
function X = Rand_NBP_Matrix(Matrix,para)
Matrices ={'NBP';
    'NBP_sequential';
    'GNBP';
    'GNBP_sequential';
    'BNBP';
    'BNBP_sequential'
    };

if nargin<1
    Matrix=Matrices{1};
end
if nargin<2
    if strcmp(Matrix,'NBP') || strcmp(Matrix,'NBP_sequential')
        para.J=10;
        para.c=0.5;
        para.gamma0=5;
        %        para.gamma0*(log(para.J+para.c)-log(para.c))
        %        pp = para.J/(para.J+para.c)
        %        para.gamma0*(log(para.J+para.c)-log(para.c))*(-1/log(1-pp)*pp/(1-pp))
    end
    if strcmp(Matrix,'GNBP') || strcmp(Matrix,'GNBP_sequential')
%         para.J=10;
%         para.c=1;
%         
%         pp=0.001:0.001:0.999;
%         [temp,dex]=min( abs(pp./(1-pp)./(log(1-para.J*log(1-pp))) - 10/12));
%         pp=pp(dex);
%         
%         para.p_i = zeros(1,10)+pp;
%         para.gamma0 = 10*(1-pp)/pp;
%         
%         qq = -log(1-para.p_i);
%         para.gamma0*(log(para.c+sum(qq))-log(para.c))
%         
%         pp=sum(qq)./(c+sum(qq));
%         (-1./log(1-pp).*pp./(1-pp))* (-1./log(1-para.p_i(1)).*para.p_i(1)./(1-para.p_i(1)))
        
        para.J=10;
        para.c=1;
        para.gamma0=4.79;
        TargetECount = 100;
        temp=gamrnd(1*ones(1,para.J),1);
        %temp = ones(1,para.J);
        pp = TargetECount/para.gamma0/para.c*temp/sum(temp);
        para.p_i = pp./(1+pp)
        qq = -log(1-para.p_i);
        para.gamma0*(log(para.c+sum(qq))-log(para.c))
        
    end
    
    if strcmp(Matrix,'BNBP') || strcmp(Matrix,'BNBP_sequential')
%         para.J=10;
%         para.c=2;
%         para.gamma0=10;
%         
%         rr=0.001:0.001:100;
%         [temp,dex]=min( abs(rr/(para.c-1)./(psi(para.c+para.J*rr)-psi(para.c))-10/12));
%         %para.gamma0*(psi(para.c+para.J*rr)-psi(para.c))
%         rr=rr(dex);
%         
%         
%         para.r = zeros(1,10)+rr;
%         para.gamma0 = 10/rr;

        para.J=10;
        para.c=2;
        para.gamma0=2;
        TargetECount = 100;
        rr=0.0001:0.0001:1000;
        %[temp,dex]=min( abs(rr/(para.c-1)*para.gamma0-TargetECount));
        [temp,dex]=min( abs(rr/(para.c-1)*12./(psi(para.c+rr)-psi(para.c))-TargetECount));
        rr=rr(dex);
        temp=gamrnd(1*ones(1,para.J),1);
        para.r = rr*temp/sum(temp);
        para.gamma0 = 12/(psi(para.c+rr)-psi(para.c))
    end
    
end



if strcmp(Matrix,'NBP') || strcmp(Matrix,'NBP_sequential')
    J = para.J;
    c = para.c;
    gamma0 = para.gamma0;
    
    if strcmp(Matrix,'NBP')
        K=poissrnd(gamma0*(log(J+c)-log(c)));
        X = zeros(J,K);
        for k=1:K
            n_k = rand_TNB(0,J/(J+c));
            prob=1/J*ones(J,1);
            X(:,k)= (multrnd_histc(n_k,prob));
            % i
        end
    end
    
    if strcmp(Matrix,'NBP_sequential')
        X = zeros(J,0);
        for j=1:J
            p_j=1/(j+c);
            Kj = poissrnd(-gamma0*(log(1-p_j)));
            for k=1:size(X,2)
                X(j,k) = nbinrnd(sum(X(:,k)),1-p_j);
            end
            for k=1:Kj
                X(j,end+1) = rand_TNB(0,p_j);
            end
        end
        K = size(X,2);
    end
end

if strcmp(Matrix,'GNBP') || strcmp(Matrix,'GNBP_sequential')
    J = para.J;
    c = para.c;
    gamma0 = para.gamma0;
    p_i = para.p_i;
    
    sumlogpi = sum(log(max(1-p_i,realmin)));
    q=-log(max(1-p_i,realmin));
    p_prime = -sumlogpi./(c-sumlogpi);
    
    if strcmp(Matrix,'GNBP')
        K = poissrnd(gamma0*(log(c-sumlogpi)-log(c)));
        
        X = zeros(J,K);
        
        for k=1:K
            l_k = rand_TNB(0,p_prime);
            l_k = multrnd_histc(l_k, (log(max(1-p_i,realmin))/sumlogpi)');
            for j=1:J
                X(j,k)=0;
                for t=1:l_k(j)
                    X(j,k)=X(j,k)+rand_TNB(0,p_i(j));
                end
            end
        end
    end
    
    if strcmp(Matrix,'GNBP_sequential')
        X = zeros(J,0);
        for j=1:J
            Kj=poissrnd(gamma0*(log(c+sum(q(1:j)))-log(c+sum(q(1:j-1)))));
            for k=1:size(X,2)
                X(j,k) = nbinrnd(sum(X(:,k)), 1-q(j)./(c+sum(q(1:j))));
            end
            for k=1:Kj
                X(j,end+1)=rand_TNB(0,q(j)./(c+sum(q(1:j))));
            end
        end
        K=size(X,2);
        for k=1:K
            for j=1:J
                temp=X(j,k);
                X(j,k)=0;
                for t=1:temp
                    X(j,k) = X(j,k)+rand_TNB(0,p_i(j));
                end
            end
        end
    end
end


if strcmp(Matrix,'BNBP') || strcmp(Matrix,'BNBP_sequential')
    J = para.J;
    c = para.c;
    gamma0 = para.gamma0;
    r = para.r;
    
    if strcmp(Matrix,'BNBP')
        
        K = poissrnd(gamma0*(psi(c+sum(r))-psi(c)));
        X = zeros(J,K);
        for k=1:K
            n_k = rand_digamma(sum(r),c,100000);
            prob=dirrnd(r(:));
            X(:,k)= (multrnd_histc(n_k,prob));
        end
    end
    
    if strcmp(Matrix,'BNBP_sequential')
        
        
        X = zeros(J,0);
        
        for j=1:J
            Kj = poissrnd(gamma0*(psi(c+sum(r(1:j)))-psi(c+sum(r(1:j-1)))));
            for k=1:size(X,2)
                pk = betarnd(sum(X(:,k)), sum(r(1:j-1))+c);
                X(j,k) = nbinrnd(r(j), 1-pk);
            end
            for k=1:Kj
                X(j,end+1) = rand_digamma(r(j),sum(r(1:j-1))+c,100000);
            end
        end
        K = size(X,2);
    end
end

