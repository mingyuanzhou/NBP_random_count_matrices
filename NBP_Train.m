function [gamma0,c,n_dot_k,r_k_ave,r_star_ave,cAve]=NBP_Train(X,Burnin,Collection,IsBinary)
% Infer the parameters for the negative binomial process 
% random count matrix using Gibbs Sampling
% Check Equation (9) of arXiv:1404.3331
if ~exist('Burnin','var')
    Burnin = 100;
end
if ~exist('Collection','var')
    Collection=100;
end
if ~exist('IsBinary','var')
    IsBinary = false;
end

c0 = 1e-3;  
d0 = 1e-3;  
e0 = 1e-3;  
f0 = 1e-3;
c=1;
gamma0=10;

[K,J]=size(X);
% J is the number of rows, as defined in the paper
% K is the number of columns, as defined in the paper

r_k_ave=zeros(K,1);
r_star_ave = 0;
gamma0_ave = 0;
n_dot_k = full(sum(X,2));
cAve=0;
n_dot_k_Ave= n_dot_k-n_dot_k;

% if IsBinary
%     Z = X;
%     [ii,jj]=find(Z);
%     iijj=Z>0;
% end

for iter=1:Burnin+Collection
    gamma0 = gamrnd(e0 + K,1/(f0 - log(max(c/(c+J),realmin))));
    r_k = gamrnd(n_dot_k, 1/(c+J));
    r_star =gamrnd(gamma0,1/(c+J));
    c = gamrnd(c0+gamma0,1./(d0+sum(r_k)+r_star));
    
    if iter>Burnin
        r_k_ave=r_k_ave+r_k;
        r_star_ave = r_star_ave+r_star;
        gamma0_ave = gamma0_ave+gamma0;
        cAve=cAve+c;
        n_dot_k_Ave =  n_dot_k_Ave +  n_dot_k;
    end
    
    %     subplot(1,3,1);plot(1:K,r_k)
    %     subplot(1,3,2);plot(mean(Z,2))
    %     %subplot(1,3,2);plot(r_i)
    %     subplot(1,3,2);imagesc(X);%colorbar
    %     drawnow
    
    %     if mod(iter,100)==0
    %         %subplot(2,1,1)
    %         %plot(1:J,p_i0,1:J,p_i,'r');
    %         %subplot(2,1,2)
    %         %imagesc(bsxfun(@rdivide,Theta,sum(Theta,1)));
    %         plot(r_k);
    %         %title([num2str(c0),'    ',num2str(c)]);
    %         drawnow
    %         pause(0.01)
    %     end
end
r_k_ave=r_k_ave/Collection;
r_star_ave = r_star_ave/Collection;
cAve=cAve/Collection;
gamma0_ave=gamma0_ave/Collection;

n_dot_k = n_dot_k_Ave/Collection;
gamma0=gamma0_ave;
c=cAve;