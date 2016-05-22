% Infer the parameters from a training count matrix and then regenerate
% a random count matrix using these parameters
% 05/2014
% Matlab code for M. Zhou, O. Madrid, and J. G. Scott, "Negative Binomial
% Process Random Count Matrices," 2014. arXiv:1404.3331
%
% Mingyuan Zhou
%

% option
%% Main code

Burnin=1000;
Collection=1500;

addpath('data/')
load TDT2.mat
Xtrain = fea';
GroundInd = gnd;

% 
%     addpath('data/')
%     load Reuters21578.mat
%     Xtrain = fea';
%     GroundInd = gnd;
% %     Xtrain = Xtrain(:,GroundInd<=30);
% %     GroundInd = GroundInd(GroundInd<=30);


% addpath('20news-bydate/')
% load train.data
% load test.data
% test(:,1)=max(train(:,1))+test(:,1);
% train_test = [train;test];
% Xtrain =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
% load train.label
% load test.label
% GroundInd = [train;test];
% 
% X = Xtrain(:,GroundInd==1);
% X = X(:,1:fix(end*0.2));

X = Xtrain(:,GroundInd==30);
X = X(sum(X,2)>0,:);

randn('state',0)
rand('state',0)

%% NBP
clear para
tic
[gamma0,c,n_k] = NBP_Train(X,Burnin,Collection);
toc
figure
xxx=hist(sum(X,2),1:1000);
xxx=xxx/sum(xxx);
prob=Logrithmic_pdf_log(1:1000,size(X,2)/(size(X,2)+c));
plot(1:20,log(xxx(1:20)+eps),1:20,prob(1:20),'r')

para.J=size(X,2);
para.c=c;
para.gamma0=gamma0;

figure;subplot(1,2,1);imagesc(log(max(order_matrix(X)',exp(-1))),[-1,3]);
if 1
    YYY_NBP=Rand_NBP_Matrix('NBP',para)';
    subplot(1,2,2);imagesc(log(max(order_matrix(YYY_NBP)',exp(-1))),[-1,3]);
else
    YYY=Rand_NBP_Matrix('NBP_sequential',para)';
    subplot(1,2,2);imagesc(log(YYY+0.1),[-1,6]);
end
subplot(1,2,1); ylabel('Documents'); xlabel('Words')
subplot(1,2,2); ylabel('Documents'); xlabel('Words')

%% GNBP
clear para
tic
[gamma0,c,p_i,r_k,r_star,l_k,c1,q] = GNBP_Train(X,Burnin,Collection);
toc
figure
xxx=hist((sum(X,2)),1:1000);
xxx=xxx/sum(xxx);

para.J=size(X,2);
para.c=c;
para.gamma0=gamma0;
para.p_i=p_i;
para.q=q;

XX = sum(Rand_NBP_Matrix('GNBP',para)',2);
for iter=1:10
    XX = [XX;sum(Rand_NBP_Matrix('GNBP',para)',2)];
end
yyy=hist(XX,1:1000);
yyy=yyy/sum(yyy);
prob=log(yyy+eps);
plot(1:20,log(xxx(1:20)+eps),1:20,prob(1:20),'r')
%plot(1:100,xxx(1:100),1:100,yyy(1:100),'r')


figure;subplot(1,2,1);imagesc(log(max(order_matrix(X)',exp(-1))),[-1,3]);
if 1
    YYY_GNBP=Rand_NBP_Matrix('GNBP',para)';
    subplot(1,2,2);imagesc(log(max(order_matrix(YYY_GNBP)',exp(-1))),[-1,3]);
else
    YYY=Rand_NBP_Matrix('GNBP_sequential',para)';
    subplot(1,2,2);imagesc(log(max(YYY,exp(-1)))',[-1,3]);
end
subplot(1,2,1); ylabel('Documents'); xlabel('Words')
subplot(1,2,2); ylabel('Documents'); xlabel('Words')

%% BNBP
clear para
tic
[gamma0,c,r_i,p_k,p_star,n_dot_k]=BNBP_Train(X,Burnin,Collection);
toc
figure
xxx=hist((sum(X(sum(X,2)>0,:),2)),1:1000);
xxx=xxx/sum(xxx);
prob=digamma_log(1:1000,sum(r_i),c);
plot(1:20,log(xxx(1:20)+eps),1:20,prob(1:20),'r')
% 
% XX = sum(Rand_NBP_Matrix('BNBP',para)',2);
% for iter=1:10
%     XX = [XX;sum(Rand_NBP_Matrix('BNBP',para)',2)];
% end
% yyy=hist(XX,1:1000);
% yyy=yyy/sum(yyy);

para.J=size(X,2);
para.c=c;
para.gamma0=gamma0;
para.r=r_i;
figure;subplot(1,2,1);imagesc(log(max(order_matrix(X)',exp(-1))),[-1,3]);
if 1
    YYY_BNBP=Rand_NBP_Matrix('BNBP',para)';
    subplot(1,2,2);imagesc(log(max(order_matrix(YYY_BNBP)',exp(-1))),[-1,3]);
else
    YYY=Rand_NBP_Matrix('BNBP_sequential',para)';
    subplot(1,2,2);imagesc(log(max(YYY,exp(-1)))',[-1,3]);
end
subplot(1,2,1); ylabel('Documents'); xlabel('Words')
subplot(1,2,2); ylabel('Documents'); xlabel('Words')


figure
subplot(2,2,1); imagesc(log(max(order_matrix(X)',exp(-1))),[-1,3]);ylabel('Documents'); xlabel('Words'); title('(a) The observed count matrix')
subplot(2,2,2); imagesc(log(max(order_matrix(YYY_NBP)',exp(-1))),[-1,3]);ylabel('Documents'); xlabel('Words'); title('(b) A simulated NBP random count matrix')
subplot(2,2,3); imagesc(log(max(order_matrix(YYY_GNBP)',exp(-1))),[-1,3]); ylabel('Documents'); xlabel('Words'); title('(c) A simulated GNBP random count matrix')
subplot(2,2,4); imagesc(log(max(order_matrix(YYY_BNBP)',exp(-1))),[-1,3]); ylabel('Documents'); xlabel('Words'); title('(d) A simulated BNBP random count matrix')

figure
subplot(2,2,1); imagesc(round_matrix(order_matrix(X)'),[0,3]);ylabel('Documents'); xlabel('Words'); title('(a) The observed count matrix')
subplot(2,2,2); imagesc(round_matrix(order_matrix(YYY_NBP)'),[0,3]);ylabel('Documents'); xlabel('Words'); title('(b) A simulated NBP random count matrix')
subplot(2,2,3); imagesc(round_matrix(order_matrix(YYY_GNBP)'),[0,3]); ylabel('Documents'); xlabel('Words'); title('(c) A simulated GNBP random count matrix')
subplot(2,2,4); imagesc(round_matrix(order_matrix(YYY_BNBP)'),[0,3]); ylabel('Documents'); xlabel('Words'); title('(d) A simulated BNBP random count matrix')
 
colormap(flipud(colormap(gray)))


set(gcf,'papersize',[30 10])
print('-dpdf','NBP_generate2.pdf')


figure
subplot(2,2,1); imagesc(log(max(order_matrix(X)',exp(-1))),[-1,3]);ylabel('Documents'); xlabel('Words'); title('(a) The observed count matrix'); set(gca,'XTick',[0:500:size(X,1),size(X,1)])
subplot(2,2,2); imagesc(log(max(order_matrix(YYY_NBP)',exp(-1))),[-1,3]);ylabel('Documents'); xlabel('Words'); title('(b) A simulated NBP random count matrix'); set(gca,'XTick',[0:500:size(YYY_NBP,1),size(YYY_NBP,1)])
subplot(2,2,3); imagesc(log(max(order_matrix(YYY_GNBP)',exp(-1))),[-1,3]); ylabel('Documents'); xlabel('Words'); title('(c) A simulated GNBP random count matrix');set(gca,'XTick',[0:500:size(YYY_GNBP,1),size(YYY_GNBP,1)])
subplot(2,2,4); imagesc(log(max(order_matrix(YYY_BNBP)',exp(-1))),[-1,3]); ylabel('Documents'); xlabel('Words'); title('(d) A simulated BNBP random count matrix');set(gca,'XTick',[0:500:size(YYY_BNBP,1),size(YYY_BNBP,1)])

set(gcf,'papersize',[30 10])
print('-dpdf','NBP_generate2.pdf')