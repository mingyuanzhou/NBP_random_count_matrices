%Demo_NBP
%
% Matlab code for M. Zhou, O. H. M. Padilla and J. G. Scott, "Priors for Random Count Matrices
% Derived from a Family of Negative Binomial Processes," arXiv:1404.3331, to appear in
% the Journal of the American Statistical Association (Theory and Methods), 2015.
%
% Copyright (c), 2015, Mingyuan Zhou
% http://mingyuanzhou.github.io/
%
% First version 12/2013
% 2nd revision 05/2014
% 3rd revision 11/2014
%
% % The K*J matrix in the code is defined as a J*K matrix in the paper
%% Demo
state=0; %0, 1, 2, 3, 4
Percentage = 10; %20, 30, 40, 50, 60, 70, 80
IsKnowKall =false; %true
option = 'NBP';
%option = 'GNBP';
%option = 'BNBP';
dataset=1; %2, 3


%% Main code

%Use S independent Markov Chains
%Run Burnin+Collection iterations for each chain
%Collect Collection samples per chain
Burnin=2499;
Collection=1;
S = 1;

%Download the example datasets and put them into the data folder
if dataset==1
    %20news-bydate data set
    %http://qwone.com/~jason/20Newsgroups/20news-bydate-matlab.tgz
    addpath('data/20news-bydate/')
    load train.data
    load test.data
    test(:,1)=max(train(:,1))+test(:,1);
    train_test = [train;test];
    Xtrain =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
    load train.label
    load test.label
    GroundInd = [train;test];
elseif dataset==2
    %Top 30 categories in TDT2
    %http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/TDT2.mat
    addpath('data/')
    load TDT2.mat
    Xtrain = fea';
    GroundInd = gnd;
else
    %CNAE-9
    %https://archive.ics.uci.edu/ml/machine-learning-databases/00233/CNAE-9.data
    addpath('data/')
    load CNAE-9.data
    data = CNAE_9;
    Xtrain = data(:,2:end)';
    GroundInd = data(:,1);
end

NumCategory=length(unique(GroundInd));

%% Create the training/testing partition for each category
rand('state',state)
[temp,dex1]=sort(rand(1,size(Xtrain,2)));
% %dex1=size(Xtrain,2):-1:1;
% dex1=1:size(Xtrain,2);
Xtrain=Xtrain(:,dex1);
GroundInd=GroundInd(dex1);

X=cell(1,NumCategory);
Kdex=cell(1,NumCategory);
dexTest=false(size(Xtrain,2),1);
for i=1:NumCategory
    dex = find(GroundInd==i);
    Len=length(dex);
    X{i} = Xtrain(:,dex(1:round(Len*Percentage/100)));
    Kdex{i} = (sum(X{i},2)>0);
    X{i} = X{i}(Kdex{i},:);
    % dexTest = [dexTest;dex((round(Len*Percentage/100)+1):length(dex))];
    dexTest(dex((round(Len*Percentage/100)+1):length(dex)))=true;
end


%% Run the code
switch option
    case 'NBP'
        InferInd_sample = zeros(NumCategory,size(Xtrain,2),S);
        for sample = 1:S
            tic
            %% Training
            n_kC = zeros(size(Xtrain,1),NumCategory);
            JC = zeros(1,NumCategory);
            cC = zeros(1,NumCategory);
            gamma0C = zeros(1,NumCategory);
            parfor i=1:NumCategory
                %i
                %fprintf('%d,',i);
                [gamma0C(i),cC(i),n_dot_k] = NBP_Train(X{i},Burnin,Collection);
                JC(i) = size(X{i},2);
                n_kC(:,i) = full(sparse(find(Kdex{i}),1,n_dot_k,size(Xtrain,1),1));
            end
            %fprintf('\n');
            %% Testing
            InferInd_sample_temp = zeros(NumCategory,nnz(dexTest));
            parfor i=1:NumCategory
                %i
                InferInd_sample_temp(i,:) = predict_NBP_Par(Xtrain(:,dexTest),n_kC(:,i),gamma0C(i),JC(i),cC(i),IsKnowKall);
                %fprintf('%d,',i);
            end
            InferInd_sample(:,dexTest,sample) = InferInd_sample_temp;
            fprintf('\nsample = %d %.0f\n',sample,toc);
        end
    case 'GNBP'
        tic
        InferInd_sample= zeros(NumCategory,size(Xtrain,2),S);
        for sample = 1:S
            %% Training
            output=cell(1,NumCategory);
            tic
            parfor i=1:NumCategory
                %i
                %fprintf('%d,',i);
                [gamma0,c,p_i,r_k,r_star,l_k,c1,q] = GNBP_Train(X{i},Burnin,Collection);
                output{i}.gamma0 = gamma0;
                output{i}.c = c;
                output{i}.p_i = p_i;
                output{i}.r_k = r_k;
                output{i}.r_star = r_star;
                output{i}.l_k = l_k;
                output{i}.q = q;
            end
            %toc
            %fprintf('\n');
            %% Testing
            
            %InferInd = zeros(NumCategory,size(Xtrain,2));
            
            LogF=LogFmatrix(max(max(Xtrain(:,dexTest))));
            
            InferInd_sample_temp = zeros(NumCategory,nnz(dexTest));
            parfor i=1:NumCategory
                %for i=1:NumCategory
                %fprintf('%d,',i);
                InferInd_sample_temp(i,:) = predict_GNBP_Par(Xtrain(:,dexTest),output{i},Kdex{i},LogF,IsKnowKall);
            end
            InferInd_sample(:,dexTest,sample)=InferInd_sample_temp;
            fprintf('\nsample = %d %.0f\n',sample,toc);
        end
    case 'BNBP'
        
        tic
        InferInd_sample = zeros(NumCategory,size(Xtrain,2),S);
        for sample = 1:S
            %% Training
            
            output=cell(1,NumCategory);
            tic
            parfor i=1:NumCategory
                %i
                %fprintf('%d,',i);
                [gamma0,c,r_i,p_k,p_star,n_dot_k]=BNBP_Train(X{i},Burnin,Collection);
                output{i}.gamma0 = gamma0;
                output{i}.c = c;
                output{i}.r_i = r_i;
                output{i}.p_k = p_k;
                output{i}.p_star = p_star;
                output{i}.n_dot_k = n_dot_k;
            end
            %toc
            %fprintf('\n');
            %% Testing
            %InferInd = zeros(NumCategory,size(Xtrain,2));
            %tic
            InferInd_sample_temp = zeros(NumCategory,nnz(dexTest));
            parfor i=1:NumCategory
                %i
                %fprintf('%d,',i);
                InferInd_sample_temp(i,:) = predict_BNBP_Par(Xtrain(:,dexTest),output{i},Kdex{i},IsKnowKall);
            end
            InferInd_sample(:,dexTest,sample)=InferInd_sample_temp;
            fprintf('sample = %d %.0f\n',sample,toc);
        end
end

% for j=1:size(Xtrain,2)
%     InferInd(:,j) = exp(InferInd(:,j) - max(InferInd(:,j)));
%     InferInd(:,j) = InferInd(:,j)/sum(InferInd(:,j));
% end




%% Accuracy using a single sample

InferInd = InferInd_sample(:,:,1);
for j=1:size(Xtrain,2)
    InferInd(:,j) = exp(InferInd(:,j) - max(InferInd(:,j)));
    InferInd(:,j) = InferInd(:,j)/sum(InferInd(:,j));
end

counts=zeros(1,5);
countalls=0;
for i=1:NumCategory
    dex = find(GroundInd==i);
    for j=1:size(Xtrain,2)
        if dexTest(j) && GroundInd(j)==i
            [temp,label]=sort(InferInd(:,j),'descend');
            for t=1:min(5,NumCategory)
                counts(t) =counts(t)+sum(label(1:t)==i);
            end
            countalls=countalls+1;
        end
    end
end
fprintf('single sample: \t \t %.4f %.4f %.4f %.4f %.4f \n ',counts/countalls)


%% Accuracy using S samples
InferInd = zeros(NumCategory,size(Xtrain,2));
for j=1:size(Xtrain,2)
    InferInd(:,j) = mean(exp(InferInd_sample(:,j,:) - max(max(InferInd_sample(:,j,:)))),3);
    InferInd(:,j) = InferInd(:,j)/sum(InferInd(:,j));
end

count=zeros(1,5);
countall=0;
for i=1:NumCategory
    dex = find(GroundInd==i);
    for j=1:size(Xtrain,2)
        if dexTest(j) && GroundInd(j)==i
            [temp,label]=sort(InferInd(:,j),'descend');
            for t=1:min(5,NumCategory)
                count(t) =count(t)+sum(label(1:t)==i);
            end
            countall=countall+1;
        end
    end
end
fprintf('%d samples:\t \t %.4f %.4f %.4f %.4f %.4f \n ',S, count/countall)


temp = InferInd;
temp(:,dex1)= InferInd;
figure;imagesc(temp)

%Vote
InferInd = zeros(NumCategory,size(Xtrain,2));
for sample = 1:S
    for j=1:size(Xtrain,2)
        [~,dex]= max(InferInd_sample(:,j,sample));
        InferInd(dex,j) = InferInd(dex,j)+1/S;
    end
end

count1=zeros(1,5);
countall1=0;
for i=1:NumCategory
    dex = find(GroundInd==i);
    for j=1:size(Xtrain,2)
        if dexTest(j) && GroundInd(j)==i
            [temp,label]=sort(InferInd(:,j),'descend');
            for t=1:min(5,NumCategory)
                count1(t) =count1(t)+sum(label(1:t)==i);
            end
            countall1=countall1+1;
        end
    end
end

fprintf('%d samples by Vote: \t %.4f %.4f %.4f %.4f %.4f \n ',S, count1/countall1)

save([option, num2str(Percentage),'_',num2str(state),'_',num2str(IsKnowKall),'_',num2str(dataset),'_samples.mat'], 'count','countall','count1','countall1','counts','countalls')