%%

clc;clear all;close all

k = 256;
urm = zeros(k);
jom = zeros(k);
jtm = zeros(k);

ss = 5;
ts = 5;

load 'C:\Users\arafat\Dropbox\CS US\cyst_d40x40_c5s40_nrf.mat'
D1 = 256.*(rf./max(rf(:)));
D1 = downsample(D1,2);

% Uniform random
for g = 1:k
    p = randperm(k);
    r = sort(p(1:floor(k/(ss*ts))));
    urm(r,g) = 1;
end

U1 = fftshift(10*log10(abs(fft2(urm))));
U1 = U1 + min(abs(U1(:)));
figure;imagesc(1-urm);colormap gray;axis square;
figure;imagesc(-(k/2)-1:k/2,-(k/2)-1:k/2,U1);colormap gray;axis square;

% 1D Jitter
for g = 1:k
    r = jitter(k,1/(ss*ts));
    jom(r,g) = 1;
end

U2 = fftshift(10*log10(abs(fft2(jom))));
U2 = U2 + min(abs(U2(:)));
figure;imagesc(1-jom);colormap gray;axis square;
figure;imagesc(-(k/2)-1:k/2,-(k/2)-1:k/2,U2);colormap gray;axis square;

% 2D Jitter
L2 = jitter2Ds(k,k,1/ts,1/ss);
jtm = jtm + L2;
U3 = fftshift(10*log10(abs(fft2(jtm))));
U3 = U3 + min(abs(U3(:)));
figure;imagesc(1-jtm);colormap gray;axis square;
figure;imagesc(-(k/2)-1:k/2,-(k/2)-1:k/2,U3);colormap gray;axis square;

