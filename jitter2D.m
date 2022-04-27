function L = jitter2D(n,m,p)
% Given integers n and m, this script subsamples the set ([n],[m]) by p in 
% a jittered manner which controls the maximum gap between missing entries.


step = round(1/p);
grid = 1:step:n;
L = zeros(n,m);
dev = floor(step/2);
for k = 1:length(grid)
    for l = 1:length(grid)
        if k ==1 || k == length(grid)
            p1 = (sign(2-k))*randperm(dev);
            p1 = p1(1);
        else
            p1 = (2*(rand<=.5) - 1)*randperm(dev);
            p1 = p1(1);
        end
        if l ==1 || l == length(grid)
            p2 = (sign(2-l))*randperm(dev);
            p2 = p2(1);
        else
            p2 = (2*(rand<=.5) - 1)*randperm(dev);
            p2 = p2(1);
        end
        L(grid(k)+p1,grid(l)+p2)=1;
    end
end

end