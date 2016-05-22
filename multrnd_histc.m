function r = multrnd_histc(n,p)
p=p(:);
edges = [0; cumsum(p,1)];
r = histc(rand(n,1)*edges(end),edges);
r(end)=[];
r=r(:);
    