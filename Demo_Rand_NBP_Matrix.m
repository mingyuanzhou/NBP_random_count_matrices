Matrices ={'NBP_row';
    'NBP_column';
    'GNBP_row';
    'GNBP_column';
    'BNBP_row';
    'BNBP_column'
    };
options={'NBP';
    'GNBP';
    'BNBP'};

count=zeros(1,100);
K = zeros(1,100);
%for i=1:100
rand('state',0)
randn('state',0)
X = Rand_NBP_Matrix('GNBP_sequential');
%X = Rand_NBP_Matrix('GNBP');
%count(i)=sum(X(:));
%K(i)=size(X,2);
%end

% subplot(1,2,1)
% 
% [x,y]=meshgrid(1:size(X,2),1:size(X,1));
% imagesc(log(X+1))
% colormap(flipud(gray))
% %colormap((gray))
% text(x(find(X(:)))-0.5,y(find(X(:))),num2str(X(find(X(:)))),'color','r') %,'HorizontalAlignment','center')

%subplot(1,2,2)


figure
rand('state',0)
randn('state',0)
X = Rand_NBP_Matrix('NBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',18)
xlabel('columns','FontSize',14)
ylabel('rows','FontSize',14)



figure
rand('state',0)
randn('state',0)
X = Rand_NBP_Matrix('GNBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',18)
xlabel('columns','FontSize',14)
ylabel('rows','FontSize',14)

figure
rand('state',0)
randn('state',0)
X = Rand_NBP_Matrix('BNBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',18)
xlabel('columns','FontSize',14)
ylabel('rows','FontSize',14)

% export_fig('NBP.pdf')
% export_fig('GNBP.pdf')
% export_fig('BNBP.pdf')


figure
for i=1:3
subplot(3,3,i)
rand('state',i)
randn('state',i)
X = Rand_NBP_Matrix('NBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',14)
xlabel('columns','FontSize',12)
ylabel('rows','FontSize',12)
title('NBP','FontSize',12)
end
for i=1:3
    subplot(3,3,3+i)
rand('state',i)
randn('state',i)
X = Rand_NBP_Matrix('GNBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',14)
xlabel('columns','FontSize',12)
ylabel('rows','FontSize',12)
title('GNBP','FontSize',12)
end
for i=1:3
      subplot(3,3,6+i)
    rand('state',i)
randn('state',i)
X = Rand_NBP_Matrix('BNBP_sequential');
[x,y]=meshgrid(1:size(X,2),1:size(X,1));
imagesc(X)
%set(gca,'dataAspectRatio',[1 1 1])
colormap(flipud(gray))
%colormap((gray))
text(x(find(X(:))),y(find(X(:))),num2str(X(find(X(:)))),'color','m','HorizontalAlignment','center','FontSize',14)
xlabel('columns','FontSize',12)
ylabel('rows','FontSize',12)
title('BNBP','FontSize',12)
end

count=0
for i=1:100
    X = Rand_NBP_Matrix('BNBP_sequential');
    count=count+sum(X(:))/100;
end
count

set(gcf,'papersize',[30 10])
print('-dpdf','NBP_Matrix_Draw.pdf')
