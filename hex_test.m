clear all;
close all;
clc;

getd = @(p)path(p,path);
getd('C:\Users\MohammadArafat\Dropbox\Elastography');

%% Data Loading
load 'org.mat'  %D1
%D1 = D1./max(D1(:));
load 'curv.mat' %recon_c
%recon_c = recon_c./max(recon_c(:));
load 'wave.mat' %recon_w
%recon_w = recon_w./max(recon_w(:));

%% Processing
env_c = envelope(recon_c);
[FX,FY] = gradient(env_c);
Gmap = abs(FX) + abs(FY); clear FX, clear FY;
Gmap = Gmap./max(Gmap(:));

recon = (1-Gmap).*recon_w + Gmap.*recon_c;

%% Evaluation
err_c = mean(abs(D1(:) - recon_c(:)).^2)
err_w = mean(abs(D1(:) - recon_w(:)).^2)
err   = mean(abs(D1(:) - recon(:)).^2)

env = (envelope(D1));
env_c = (envelope(recon_c));
env_w = (envelope(abs(recon_w)));
env_com = (envelope(recon));

option = 'mu';
fac = 255;
y = LogComp(env,option,fac);
y = 1 + y./max(y(:));
y1 = LogComp(env_c,option,fac);
y1 = 1 + y1./max(y1(:));
y2 = LogComp(env_w,option,fac);
y2 = 1 + y2./max(y2(:));
y3 = LogComp(env_com,option,fac);
y3 = 1 + y3./max(y3(:));


MSE_c = sum(sum((y-y1).^2))/size(y1,1)/size(y1,2)
MSE_w = sum(sum((y-y2).^2))/size(y2,1)/size(y2,2)
MSE = sum(sum((y-y3).^2))/size(y3,1)/size(y3,2)


[ssim1, ssim_map] = ssim(y, y1, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM_c = ssim1
[ssim2, ssim_map] = ssim(y, y2, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM_w = ssim2
[ssim3, ssim_map] = ssim(y, y3, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM = ssim3

figure; imagesc(y);colormap gray;axis off;%title('Actual')
figure; imagesc(y1);colormap gray;axis off;%title('Curvelet')
figure; imagesc(y2);colormap gray;axis off;%title('Waveatom')
figure; imagesc(y3);colormap gray;axis off;%title('Combined')