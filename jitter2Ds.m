function L = jitter2Ds(n,m,p,q)


step1 = round(1/p);
step2 = round(1/q);
grid1 = 1:step1:n;
grid2 = 1:step2:m;
L = zeros(n,m);
dev1 = floor(step1/2);
dev2 = floor(step2/2);

for k = 1:length(grid1)
    for l = 1:length(grid2)
        if k == 1 || k == length(grid1)
            p1 = (sign(2-k))*randperm(dev1);
            p1 = p1(1);
        else
            p1 = (2*(rand<=.5) - 1)*randperm(dev1);
            p1 = p1(1);
        end
        if l == 1 || l == length(grid2)
            p2 = (sign(2-l))*randperm(dev2);
            p2 = p2(1);
        else
            p2 = (2*(rand<=.5) - 1)*randperm(dev2);
            p2 = p2(1);
        end
        L(grid1(k)+p1,grid2(l)+p2)=1;
    end
end

end