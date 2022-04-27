%%
close all;
clear all;
clc;

getd = @(p)path(p,path);
getd('C:\Users\arafat\Dropbox\Elastography\spot-1.2');

% Original data D
%load 'Z:\Field_II_new\cyst_phantom_RF.mat';
%load 'C:\Users\arafat\Dropbox\Elastography\bmode_spine.mat';
%D = 100.*image_data./max(image_data(:));
D = bmode_spine;
[nr,ns] = size(D);

% Original data spectrum
F = opDFT2(nr,ns);
FD = fftshift(reshape(F*D(:),nr,ns));

% Define undersampling percentage
ss = 12;
ts = 8;
perc = 1/ss;
perct = 1/ts;

% Regular sampling
Drr = zeros(size(D));
Drr(1:ss:nr,1:ts:ns) = D(1:ss:nr,1:ts:ns);
FDrr = fftshift(reshape(F*Drr(:),nr,ns));
Rx = opRestriction(nr,1:ss:nr);
Ry = opRestriction(ns,1:ts:ns);
Rreg = opKron(Rx,Ry);

% Uniform random sampling along scan-line
p = randperm(nr);
r = sort(p(1:floor(nr/ss)));
Dir = zeros(size(D));
Dir(r,1:ts:ns) = D(r,1:ts:ns);
FDir = fftshift(reshape(F*Dir(:),nr,ns));
Rx = opRestriction(nr,r);
Ry = opRestriction(ns,1:ts:ns);
Runf = opKron(Rx,Ry);

% Uniform random sampling along transducer
p = randperm(ns);
r = sort(p(1:floor(ns/ts)));
Dir = zeros(size(D));
Dir(1:ss:nr,r) = D(1:ss:nr,r);
FDirt = fftshift(reshape(F*Dir(:),nr,ns));
Rx = opRestriction(nr,1:ss:nr);
Ry = opRestriction(ns,r);
Runt = opKron(Rx,Ry);

% Jitter random sampling along scan-line
r = jitter(nr,perc);
Djr = zeros(size(D));
Djr(r,1:ts:ns) = D(r,1:ts:ns);
FDjr = fftshift(reshape(F*Djr(:),nr,ns));
Rx = opRestriction(nr,r);
Ry = opRestriction(ns,1:ts:ns);
Rjit1 = opKron(Rx,Ry);

% Jitter random sampling along trasducer
r = jitter(ns,perct);
Djr = zeros(size(D));
Djr(1:ss:nr,r) = D(1:ss:nr,r);
FDjrt = fftshift(reshape(F*Djr(:),nr,ns));
Rx = opRestriction(nr,1:ss:nr);
Ry = opRestriction(ns,r);
Rjit2 = opKron(Rx,Ry);

% Jitter random sampling along scan-line and transducer (separable)
ax1 = jitter(nr,perc);
L = zeros(nr,ns);
for k = 1:length(ax1)
    ax2 = jitter(ns,perct);
    L(ax1(k),ax2) = 1;
end
Rjin1 = opMask(ns*nr,find(L(:))); % Restriction Operator
FDnjj = fftshift(reshape(F*Rjin1'*(Rjin1*D(:)),nr,ns)); % Apply Restriction + Fourier and reshape

% Full Jitter random sampling along scan-line and transducer (nonseparable)
L2 = jitter2Ds(ns,nr,perc,perct);
Rjin2 = opMask(ns*nr,find(L2(:))); % Restriction Operator
FDnjnj = fftshift(reshape(F*Rjin2'*(Rjin2*D(:)),nr,ns)); % Apply Restriction + Fourier and reshape

% Plot the respective spectrums

figure;
imagesc(abs(FDnjnj));title('Full Jitter random sampling along scan-line and transducer (nonseparable)');
figure;
imagesc(abs(FDnjj));title('Jitter random sampling along scan-line and transducer (separable)');
figure;
imagesc(abs(FDjr));title('Jitter random sampling along scan-line');
figure;
imagesc(abs(FDir));title('Uniform random sampling along scan-line');
figure;
imagesc(abs(FDirt));title('Uniform random sampling along transducer');
figure;
imagesc(abs(FDjrt));title('Jitter random sampling along transducer');
figure;
imagesc(abs(FDrr));title('Regular sampling');
figure;
imagesc(abs(FD));title('Original');

