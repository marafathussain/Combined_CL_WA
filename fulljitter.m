function D = fulljitter(A,mode)
if nargin == 1
    mode = 1;
end
y = mode;
D1 = reshape(A,354,354);
load Fjitter
xir = Fjitter(1,:);
yir = Fjitter(2,:);

D = zeros(size(D1));
for k = 1:length(xir)
    D(xir(k),yir(k)) = D1(xir(k),yir(k));
end

D = D(:);

end