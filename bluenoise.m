clear all;
close all;
clc;


%% 1D
N=266;
n=38;

% Su=ignal
y=zeros(N,1);

%  No large gap random sampling
u=randi(3,38,1);
indr=((1:N/n:N)'+randi(N/n-1,38,1))';

yr=y;
yr(indr)=1; % Jitter Sampling

% Large gap randon sampling

indi=randperm(100);
yi=y;
yi(indi(1:floor(n/2)+1))=1;
yi(indi(1:floor(n/2))+150)=1; %uniform random

% Fourier transform and frequency range
F=opDFT(N);
freq=linspace(-pi,pi,N)';

% Fourier transform and DC component removal
fr=F*yr;fr(1)=0;
fi=F*yi;fi(1)=0;

% Sampling plot
plot(yr,'-or');hold on;plot(yi,'-*b');hold off;axis([0 250 0 2]);
legend('regular random sampling','rarge gap sampling');
% Centered plot
figure;
plot(freq,abs(fftshift(fr)),'r');hold on;
plot(freq,abs(fftshift(fi)),'b');hold off;
axis([-pi pi 0 1]);
ylabel('spectral magnitude');
xlabel('frequency [rad]');
legend('regular sampling','large gap sampling');

%% 2D


%% 
% In the paper it's made by kroenecker product, 
% I just saw it and it's propably easier and samrter for large matrix, I'll look into it


N=266;
n=38;

% Signal
y=zeros(N,N);

% Regular

yr=y;

ind =7:7:N-7;
yr(ind,ind)=1;

figure;
subplot(3,2,1);
imagesc(yr,[0,1]);colormap(gray);
xlabel('Receiver No');ylabel('Source No');

% Uniform random along one axis, regular along second axis

yrr=y;
randind=randperm(N);
indr=randind(1:N/7);

yrr(ind,indr)=1;

subplot(3,2,2);
imagesc(yrr,[0,1]);
xlabel('Receiver No');ylabel('Source No');

% jitteres along one axis, regular along second axis

yjr=y;

indj=min(max((1:7:N)'+randi(7,N/7,1),1),N)';

yjr(ind,indj)=1;

subplot(3,2,3);
imagesc(yjr,[0,1]);
xlabel('Receiver No');ylabel('Source No');
% 2D jittered

yjj=y;

indj1=min(max((1:7:N)'+randi(7,N/7,1),1),N)';
indj2=min(max((1:7:N)'+randi(7,N/7,1),1),N)';
yjj(indj1,indj2)=1;

subplot(3,2,4);
imagesc(yjj,[0,1]);
xlabel('Receiver No');ylabel('Source No');
% 2D jittered, different jittering along each axis


yj2=y;
indj1=min(max((1:7:N)'+randi(7,N/7,1),1),N)';

for i=1:length(indj1);
    indj2=min(max((1:7:N)'+randi(7,N/7,1),1),N)';
    yj2(indj2,indj1(i))=1;
end

subplot(3,2,5);
imagesc(yj2,[0,1]);
xlabel('Receiver No');ylabel('Source No');
% fully 2D jittered

yfj=y;

for i=1:length(ind)
    for j=1:length(ind)
        k=ind(i)+randi([-7,7]);
        l=ind(j)+randi([-7,7]);
        yfj(min(max(k,1),N),min(max(l,1),N))=1;
    end
end

subplot(3,2,6);
imagesc(yfj,[0,1]);
xlabel('Receiver No');ylabel('Source No');

%% Fourier of 2D

F=opDFT2(N,N,1);

fr=fftshift(reshape(abs(F*yr(:)),N,N));
frr=fftshift(reshape(abs(F*yrr(:)),N,N));
fjr=fftshift(reshape(abs(F*yjr(:)),N,N));
fjj=fftshift(reshape(abs(F*yjj(:)),N,N));
fj2=fftshift(reshape(abs(F*yj2(:)),N,N));
ffj=fftshift(reshape(abs(F*yfj(:)),N,N));

figure;
subplot(3,2,1);
imagesc(fr,[0,1]);
subplot(3,2,2);
imagesc(frr,[0,1]);
subplot(3,2,3);
imagesc(fjr,[0,1]);
subplot(3,2,4);
imagesc(fjj,[0,1]);
subplot(3,2,5);
imagesc(fj2,[0,1]);
subplot(3,2,6);
imagesc(ffj,[0,1]);