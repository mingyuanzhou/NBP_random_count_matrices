%Multinomial naive Bayes classifier with Laplace smoothing
%
% This is the demo code for 20newsgroups using the defaul by-date
% traning/testing partition.

addpath('data/20news-bydate/')
load train.data
load test.data
test(:,1)=max(train(:,1))+test(:,1);
train_test = [train;test];
Xtrain =sparse(train_test(:,2),train_test(:,1),train_test(:,3));
load train.label
load test.label
GroundInd = [train;test];

NumCategory=length(unique(GroundInd));
rand('state',state)
randn('state',state)
[temp,dex1]=sort(1:size(Xtrain,2));
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
    X{i} = Xtrain(:,dex(dex<=length(train)));
    Kdex{i} = (sum(X{i},2)>0);
    X{i} = X{i}(Kdex{i},:);
    % dexTest = [dexTest;dex((round(Len*Percentage/100)+1):length(dex))];
    dexTest(dex(dex>length(train)))=true;
end


alpha=1;


InferInd = zeros(NumCategory,size(Xtrain,2));
%InferInd = InferInd0;


for i=1:NumCategory
    %i
    
    prob=zeros(size(Xtrain,1),1);
    prob(Kdex{i}) = sum(X{i},2);
    prob = prob+alpha;
    
    logprob = log(prob/sum(prob));
    
    parfor j=1:size(Xtrain,2)
        
        if dexTest(j)
            InferInd(i,j)= sum(Xtrain(:,j).*logprob); %/sum(Xtrain(:,j));
        end
    end

end


InferInd0=InferInd;

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


temp = InferInd;
temp(:,dex1)= InferInd;



