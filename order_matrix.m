function [Y,dex2] = order_matrix(X)
dex=1:size(X,1);
dex2=[];
for j=1:size(X,2)
    dex1=find(X(dex,j));
    dex2=[dex2,dex(dex1)];
    dex(dex1)=[];
end
Y=X(dex2,:);