function grid = jitter(n,p)
% Given integer n, this script subsamples the set [n] by p in a jittered
% manner which controls the maximum gap between missing entries.


step = round(1/p);
grid = 1:step:n;
dev = floor(step/2);
for k = 1:length(grid)
    if k == 1 || k == length(grid)
    p = (sign(2-k))*randperm(dev);
    grid(k) = grid(k)+p(1);
    else
    p = (2*(rand<=.5) - 1)*randperm(dev);
    grid(k) = grid(k)+p(1);
    end
end

end