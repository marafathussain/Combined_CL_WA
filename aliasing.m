close all;
clear all;
clc;

getd = @(p)path(p,path);
getd('C:\Users\arafat\Dropbox\Elastography\spot-1.2');

% Fully Sampled data D
load rf1.mat;
D = rf1;
[nr,ns] = size(D);

% Fully Sampled data spectrum
F = opDFT2(nr,ns);
FD = fftshift(reshape(F*D(:),nr,ns));

% Regular sampling
Drr = zeros(size(D));
Drr(1:2:nr,1:2:ns) = D(1:2:nr,1:2:ns);
FDrr = fftshift(reshape(F*Drr(:),nr,ns));
Rx = opRestriction(nr,1:2:nr);
Ry = opRestriction(ns,1:2:ns);
Rreg = opKron(Rx,Ry);

% Uniform random sampling along receiver
p = randperm(nr);
r = sort(p(1:nr/2));
Dir = zeros(size(D));
Dir(1:2:nr,r) = D(1:2:nr,r);
FDir = fftshift(reshape(F*Dir(:),nr,ns));
Rx = opRestriction(nr,1:2:nr);
Ry = opRestriction(ns,r);
Runf = opKron(Rx,Ry);

% Jitter random sampling along receiver
r = jitter(nr);
Djr = zeros(size(D));
Djr(1:2:nr,r) = D(1:2:nr,r);
FDjr = fftshift(reshape(F*Djr(:),nr,ns));
Rx = opRestriction(nr,1:2:nr);
Ry = opRestriction(ns,r);
Rjit1 = opKron(Rx,Ry);

% Jitter random sampling along receiver and sources
r = jitter(nr);
r2 = jitter(nr);
Djj = zeros(size(D));
Djj(r2,r) = D(r2,r);
FDjj = fftshift(reshape(F*Djj(:),nr,ns));
Rx = opRestriction(nr,r);
Ry = opRestriction(ns,r2);
Rjit2 = opKron(Rx,Ry);

% Jitter random sampling along receiver and sources (nonseparable)
r = jitter(nr);
l = length(r);
R = zeros(l+1,l);
R(1,:) = r;
Dnjj = zeros(size(D));
%Rjin1 = 0;
for k = 1:l
r2 = jitter(nr);
R(k+1,:) = r2';
Dnjj(r2,r(k)) = D(r2,r(k));
end
FDnjj = fftshift(reshape(F*Dnjj(:),nr,ns));
Rjin1 = opFunction(177*177,354*354,@semijitter); 

% Full Jitter random sampling along receiver and sources (nonseparable)
[xr,yr] = ndgrid(1:2:nr,1:2:ns);
dx = (2*(rand(size(xr))<=.5) -1).*(rand(size(xr))<=.5);
dy = (2*(rand(size(yr))<=.5) -1).*(rand(size(yr))<=.5);
dx(1,:) = 0;
dx(end,:) = 0;
dx(:,1) = 0;
dx(:,end) = 0;
dy(1,:) = 0;
dy(end,:) = 0;
dy(:,1) = 0;
dy(:,end) = 0;

xir = xr(:) + dx(:); % Our random grid
yir = yr(:) + dy(:);

Dnjnj = zeros(size(D));
for k = 1:length(xir)
    Dnjnj(xir(k),yir(k)) = D(xir(k),yir(k));
end
FDnjnj = fftshift(reshape(F*Dnjnj(:),nr,ns));
Rjin2 = opFunction(354*354,354*354,@fulljitter);

%Plot the respective spectrums

figure;
imagesc(abs(FDnjnj));colormap('bone');title('Full Jitter random sampling along receiver and sources (nonseparable)');
figure;
imagesc(abs(FDnjj));colormap('bone');title('Jitter random sampling along receiver and sources (nonseparable)');
figure;
imagesc(abs(FDjj));colormap('bone');title('Jitter random sampling along receiver and sources');
figure;
imagesc(abs(FDjr));colormap('bone');title('Jitter random sampling along receiver');
figure;
imagesc(abs(FDir));colormap('bone');title('Uniform random sampling along receiver');
figure;
imagesc(abs(FDrr));colormap('bone');title('Regular sampling');
figure;
imagesc(abs(FD));colormap('bone');title('Full Sampling');

