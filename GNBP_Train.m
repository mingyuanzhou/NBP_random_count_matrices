function [gamma0,c,p_i,r_k,r_star,L_k,c_ave,q_ave]=GNBP_Train(X,Burnin,Collection,IsBinary)
% Infer the parameters for the gamma-negative binomial process 
% random count matrix using Gibbs Sampling
% Check Equation (30) of arXiv:1404.3331
if ~exist('Burnin','var')
    Burnin = 100;
end
if ~exist('Collection','var')
    Collection=100;
end
if ~exist('IsBinary','var')
    IsBinary = false;
end


[K,J]=size(X);
% J is the number of rows, as defined in the paper
% K is the number of columns, as defined in the paper
r_star_ave = 0;
r_k_ave=zeros(K,1);
p_i=ones(1,J)*0.5;
c=1;
sumlogpi = sum(log(max(1-p_i,realmin)));
p_prime = -sumlogpi./(c-sumlogpi);
a0=1e-3;    b0=1e-3;    e0=1e-3;    f0=1e-3; c0=1e-3; d0=1e-3;
r_k=ones(K,1);
l_k_ave = zeros(K,1);
c_ave=0;
q_ave=0;
p_i_ave=p_i-p_i;
gamma0_ave=0;

% if IsBinary
%     Z = X;
%     [ii,jj]=find(X);
%     iijj=Z>0;
% end

XT = sparse(X');
for iter=1:Burnin+Collection
    %iter
    
    gamma0 = gamrnd(e0 + K,1/(f0 - log(max(1-p_prime,realmin))));
    L_k = CRT_sum_mex_matrix(XT,r_k')';
    %L_k=zeros(K,1);
    %parfor k=1:K
    %    L_k(k) = CRT_sum_mex((full(X(k,:))),r_k(k));
    %end
    r_k = gamrnd(L_k, 1./(-sumlogpi+ c));
    r_star = gamrnd(gamma0, 1./(-sumlogpi+ c));
    c = gamrnd(c0+gamma0,1/(d0+sum(r_k)+r_star));
    p_i = betarnd(a0 + sum(X,1),b0+sum(r_k)+r_star);
    %p_i = betarnd(a0 + sum(sum(X,1)),b0+J*(sum(r_k)+r_star))*ones(1,J);
    sumlogpi = sum(log(max(1-p_i,realmin)));
    p_prime = -sumlogpi./(c-sumlogpi);
       
     if iter>Burnin  
        r_k_ave=r_k_ave+r_k;
        r_star_ave = r_star_ave+r_star;
        l_k_ave = l_k_ave + L_k;
        c_ave=c_ave+c;
        q_ave=q_ave-sumlogpi;
        p_i_ave=p_i_ave+p_i;
        gamma0_ave=gamma0_ave+gamma0;
     end
     
            
      %  subplot(1,4,1);plot(1:K,r_k)
%        subplot(1,4,2);plot(mean(Z,2))
      %  subplot(1,4,3);plot(p_i)
     %   subplot(1,4,4);imagesc(X);%colorbar
    %    drawnow
     
    if 0 %mod(iter,100)==0
        subplot(2,2,1)
        plot(1:J,p_i,'r');
        subplot(2,2,2)
        plot(r_k)
       % subplot(2,2,3)
       % imagesc(bsxfun(@rdivide,Theta(1:100,:),sum(Theta(1:100,:),1)));
       % title([num2str(c0),'    ',num2str(c)]);
%         subplot(2,2,4)
%         sumlogpi_new=0;
%         lambda=KK;
%         J0=size(X0,2);
%         for j=1:J0
%             if j>J
%                 pj=betarnd(1e-2+sum(X0(:,j)), 1e-2+sum(r_k)+r_star);
%             else
%                 pj=p_i(j);
%             end
%             sumlogpi_new = sumlogpi_new + log(1-pj);
%             lambda(j) = gamma0*(log(c-sumlogpi_new)-log(c-0*sumlogpi));
%         end
        %lambda=gamma0*(log(c-sumlogpi_new)-log(c-sumlogpi));
        %      plot(size(X0,1)-K+round((-100:100)*lambda/100),poisspdf(size(X0,1)-K+round((-100:100)*lambda/100),lambda));hold on
        %     stem(size(X0,1)-K,0.01);hold off
        %plot(1:J0,KK,1:J0,lambda,'r')
        drawnow
        pause(0.01)
    end
end
r_k_ave = r_k_ave/Collection;
r_star_ave = r_star_ave/Collection;
l_k_ave = l_k_ave/Collection;
c_ave=c_ave/Collection;
q_ave=q_ave/Collection;
gamma0_ave=gamma0_ave/Collection;
p_i_ave=p_i_ave/Collection;

gamma0 = gamma0_ave;
c = c_ave;
p_i = p_i_ave;
r_k = r_k_ave;
r_star = r_star_ave;
L_k = l_k_ave;