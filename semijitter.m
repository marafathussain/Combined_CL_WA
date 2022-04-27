function D = semijitter(A,mode)
if nargin ==1
    mode =1;
end
y=mode;

D1 = reshape(A,354,354);
load Sjitter
r = Sjitter(1,:);

D = zeros(size(D1));
for k = 1:177
    r2 = Sjitter(k+1,:);
    D(r2,r(k)) = D1(r2,r(k));
end

D = D(:);

end