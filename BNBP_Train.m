function [gamma0,c,r_i,p_k,p_star,n_dot_k]=BNBP_Train(X,Burnin,Collection,IsBinary)
% Infer the parameters for the beta-negative binomial process 
% random count matrix using Gibbs Sampling
% Check Equation (34) of arXiv:1404.3331
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
p_star_ave = 0;
p_k_ave=zeros(K,1);
c_ave = 0;
gamma0_ave=0;
r_i_ave=zeros(1,J);
n_dot_k = sum(X,2);
n_dot_k_Ave =  zeros(K,1);

c=1;
%sumlogpi = sum(log(max(1-p_i,realmin)));
%p_prime = -sumlogpi./(c-sumlogpi);
a0=1e-3;    b0=1e-3;    e0=1e-3;    f0=1e-3;
%r_k=ones(K,1);
%Theta =zeros(K,J);
%ThetaSum=zeros(K,J);
%l_k_ave = zeros(K,1);
%c_ave=0;
%q_ave=0;

c0=1e-0;
d0=1e-0;


gamma0=1;
c=1;
c0=c;
r_i=ones(1,J);
p_star=0;
count=0;



for iter=1:Burnin+Collection
    
    %gamma0 = gamrnd(1e-2+K,1./(1e-2+c*(psi(c+sum(r_i))-psi(c))));
    gamma0 = gamrnd(e0+K,1./(f0+psi(c+sum(r_i))-psi(c)));
    p_k = betarnd(full(sum(X,2)),c+sum(r_i));
    %if iter>10
    % p_star = log_betaProcess_rnd(1,gamma0*c,c+sum(r_i));
    p_star = logBeta_rnd(1,gamma0,c+sum(r_i));
    %p_star = gamma0*c*psi(1,c+sum(r_i));
    %end
    L_i = CRT_sum_mex_matrix(sparse(X),r_i);
    sumlogpk = sum(log(max(1-p_k,realmin)));
    r_i = gamrnd(a0 + L_i, 1./(-sumlogpk + p_star + b0));
    %r_i = gamrnd(a0 + sum(L_i), 1./((-sumlogpk + p_star)*J + b0))*ones(1,J);
    % c = gamrnd(1e-2+gamma0,1./(1e-2 -sumlogpk + p_star));
    
    if 1
        %          cnew = gamrnd(1e-0+ K,1./(1e-0+ gamma0*(psi(sum(r_i)+c)-psi(c)) ));
        %
        %          temp = log(gampdf(c, 1e-0+ K, 1./(1e-0+ gamma0*(psi(sum(r_i)+cnew)-psi(cnew)) )));
        %          temp = temp - log(gampdf(cnew, 1e-0+ K, 1./(1e-0+ gamma0*(psi(sum(r_i)+c)-psi(c)))));
        %        %  temp = temp + (-gamma0*cnew*(psi(sum(r_i)+cnew)-psi(cnew)) +cnew*sumlogpk + nnz(sum(ZSDS,2)>0)*log(cnew) );
        %        %  temp = temp - (-gamma0*c*(psi(sum(r_i)+c)-psi(c)) +c*sumlogpk + nnz(sum(ZSDS,2)>0)*log(c));
        %
        %
        %         temp = temp + (-gamma0*cnew*(psi(sum(r_i)+cnew)-psi(cnew)) +...
        %              K*log(cnew+sum(r_i)) - sum(log(cnew+sum(r_i)+sum(X,2)))  + K*log(cnew) );
        %          temp = temp - (-gamma0*c*(psi(sum(r_i)+c)-psi(c)) +...
        %              K*log(c+sum(r_i)) - sum(log(c+sum(r_i)+sum(X,2)))  + K*log(c));
        
        sumlogpk = -sum(p_k);
        cnew = gamrnd(c0+gamma0,1./(d0 -sumlogpk + 1*p_star));% -sumlogpk + p_star));
        %cnew = max(c+randn(1)*0.1,0);
        
        
        
        temp = (-gamma0*(psi(sum(r_i)+cnew)-psi(cnew)) +...
            K*gammaln(cnew+sum(r_i)) - sum(gammaln(cnew+sum(r_i)+sum(X,2)))  + 0*K*log(cnew) ) ;
        %temp = temp+ log(gampdf(cnew,1e-0,1./1e-0));
        temp = temp + (c0-1)*log(cnew)-d0*cnew;
        temp = temp - (-gamma0*(psi(sum(r_i)+c)-psi(c)) +...
            K*gammaln(c+sum(r_i)) - sum(gammaln(c+sum(r_i)+sum(X,2)))  + 0*K*log(c));
        %temp = temp- log(gampdf(c,1e-0,1./1e-0));
        temp = temp -( (c0-1)*log(c)-d0*c );
        
        %temp = temp + log(gampdf(c, 1e-0+gamma0,1./(1e-0 -sumlogpk + 1*p_star)));
        temp = temp + (c0+gamma0-1)*log(c) - (d0 -sumlogpk + 1*p_star)*c;
        %temp = temp - log(gampdf(cnew, 1e-0+gamma0,1./(1e-0 -sumlogpk + 1*p_star)));
        temp = temp - ((c0+gamma0-1)*log(cnew) - (d0 -sumlogpk + 1*p_star)*cnew);
        %             prob = exp(logprob - max(logprob));
        %             cdf = cumsum(prob);
        %             if iter>0.45
        %                 [temp,k]=max(prob);
        %             else
        %                 k = sum(rand(1)*cdf(end)>cdf)+1;
        %             end
        
        %temp
        if  exp(temp)>rand(1)
            c=cnew;
            count=count+1;
            
            %else
            %    iter
        end
        %c=cnew;
        % [count/iter,c,cnew]
    end
    
    % Theta = gamrnd(X + ones(K,1)*r_i, p_k*ones(1,J));
    

    
    
    
    if iter>Burnin
        
        p_k_ave=p_k_ave+p_k;
        p_star_ave = p_star_ave+p_star;
        gamma0_ave = gamma0_ave + gamma0;
        r_i_ave = r_i_ave + r_i;
        c_ave = c_ave +c;
        % l_k_ave = l_k_ave + L_k;
        %c_ave=c_ave+c;
        %q_ave=q_ave-sumlogpi;
        %   ThetaSum=ThetaSum+Theta;
        %   Theta=ThetaSum/(iter-Burnin);
        n_dot_k_Ave =  n_dot_k_Ave +  n_dot_k;
    end
    
%     subplot(1,3,1);plot(1:K,p_k,'r',1:K,mean(Z,2))
%     %subplot(1,3,2);plot(mean(Z,2))
%     subplot(1,3,2);plot(r_i)
%     subplot(1,3,3);imagesc(X);%colorbar
%     drawnow
    
    if 0 %mod(iter,100)==0
        subplot(2,2,1)
        plot(1:J,r_i,'r');
        title([num2str(c0),'    ',num2str(c)]);
%         subplot(2,2,3)
%         if iter<=Burnin
%             imagesc(bsxfun(@rdivide,Theta,sum(Theta,1)));
%         else
%             imagesc(bsxfun(@rdivide,ThetaSum,sum(ThetaSum,1)));
%         end
        
%         subplot(2,2,4)
%         sumr=0;
%         lambda=KK;
%         J0=size(X0,2);
%         r_i_temp = mean(r_i)*ones(1,J0-J);
%         for iiii=1:10
%             L_i_temp = CRT_sum_mex_matrix(sparse(X0(:,J+1:J0)),r_i_temp);
%             r_i_temp = gamrnd(1e-6 + L_i_temp, 1./(-sumlogpk + p_star + 1e-6));
%         end
%         
%         for j=1:J0
%             if j>J
%                 rj=r_i_temp(j-J);
%             else
%                 rj=r_i(j);
%             end
%             sumr = sumr + rj;
%             lambda(j) = gamma0*(psi(c+sumr)-psi(c));
%         end
%         %lambda=gamma0*(log(c-sumlogpi_new)-log(c-sumlogpi));
%         %      plot(size(X0,1)-K+round((-100:100)*lambda/100),poisspdf(size(X0,1)-K+round((-100:100)*lambda/100),lambda));hold on
%         %     stem(size(X0,1)-K,0.01);hold off
%         plot(1:J0,KK,1:J0,lambda,'r')
        
        drawnow
    end
    %pause(0.05)
end

p_k_ave = p_k_ave/Collection;
p_star_ave = p_star_ave/Collection;
gamma0_ave = gamma0_ave/Collection;
c_ave=c_ave/Collection;
r_i_ave = r_i_ave/Collection;

n_dot_k=n_dot_k_Ave/Collection;
gamma0 = gamma0_ave;
c = c_ave;
r_i = r_i_ave;
p_k = p_k_ave;
p_star = p_star_ave;

% p_k_ave = p_k;
% p_star_ave = p_star;
% gamma0_ave = gamma0;
% c_ave=c;
% r_i_ave = r_i;

% l_k_ave = l_k_ave/Collection;
% c_ave=c_ave/Collection;
% q_ave=q_ave/Collection;