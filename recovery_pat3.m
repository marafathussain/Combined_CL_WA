%%
% Load Data
close all;
clear all;
clc;

getd = @(p)path(p,path);
getd('C:\Users\arafat\Dropbox\WA\spot-1.2.wa');
getd('C:\Users\arafat\Dropbox\WA\spgl1-1.8.wa');
getd('C:\Users\arafat\Dropbox\WA\WaveAtom-1.1.1');

% load 'C:\Users\arafat\Dropbox\CS US\cyst_d40x40_c5s40_nrf.mat'
% D1 = 256.*(rf./max(rf(:)));
% D1 = downsample(D1,2);
load 'C:\Users\arafat\Dropbox\Elastography\AM2D\rf011.mat'
rf = RfDataFilt(1:2304,251:379);clear RfDataFilt;
rf = downsample(rf,2);
D1 = 256.*(rf./max(rf(:)));

recon_ur = zeros(size(D1));
recon_jo = recon_ur;
recon_jt = recon_ur;
[xx,yy] = size(D1);
uni_err = 0;
jit_err_1d = 0;
jit_err_2d = 0;

% Define undersampling percentage
ss = 5;
ts = 2;
perc = 1/ss;
perct = 1/ts;
k = 128;

for i = 1:k:xx-k+1
    for j = 1:k:yy-k+1
        
        urm = zeros(k);
        jom = zeros(k);

        D = D1(i:i+k-1,j:j+k-1);
        [nr,ns] = size(D);

        % Uniform random sampling along scan-line
        for g = 1:k
            p = randperm(nr);
            r = sort(p(1:floor(nr/(ss*ts))));
            urm(r,g) = 1;
        end
        Runf = opMask(ns*nr,find(urm(:))); % Restriction Operator

        % Jitter random sampling along trasducer
        for g = 1:k
            r = jitter(nr,1/(ss*ts));
            jom(r,g) = 1;
        end
        Rjit1 = opMask(ns*nr,find(jom(:))); % Restriction Operator

        % Full Jitter random sampling along scan-line and transducer (nonseparable)
        L2 = jitter2Ds(ns,nr,perct,perc);
        Rjin2 = opMask(ns*nr,find(L2(:))); % Restriction Operator

        
        % Instruct spgl1 recovery 
        C = opWaveatom(nr,ns); %Sparse Transform
        options = spgSetParms('optTol', 5e-4, 'bpTol', 5e-4, ...
                              'iterations', 100, 'verbosity', 0);


        % Uniform random sampling along scan-line
        RD = Runf*D(:);
        b = RD(:); % Observed subsampled data
        A = Runf*C'; % Measurement Operator
        xunf = spg_bp(A,b,options); %Use basis pursuit to recover our solution
        dunf = C'*xunf; %Transform recovery back to time-space domain
        Dunf = reshape(dunf,nr,ns); 
        recon_ur(i:i+k-1,j:j+k-1) = Dunf;        
        Ddiffunf = D - Dunf; %Our Error
        SNRunf = mean(abs(Ddiffunf(:)).^2); %Signal to Noise ratio of recovery
        uni_err = uni_err + SNRunf;  

        % Jitter random sampling along scan-line
        RD = Rjit1*D(:);
        b = RD(:); % Observed subsampled data
        A = Rjit1*C'; % Measurement Operator  
        xjit1 = spg_bp(A,b,options); %Use basis pursuit to recover our solution
        djit1 = C'*xjit1; %Transform recovery back to time-space domain
        Djit1 = reshape(djit1,nr,ns); 
        recon_jo(i:i+k-1,j:j+k-1) = Djit1;        
        Ddiffjit1 = D - Djit1; %Our Error
        SNRjit1 = mean(abs(Ddiffjit1(:)).^2); %Signal to Noise ratio of recovery
        jit_err_1d = jit_err_1d + SNRjit1; 

        % Full Jitter random sampling along receiver and sources (nonseparable)
        RD = Rjin2*D(:);
        b = RD(:); % Observed subsampled data
        A = Rjin2*C'; % Measurement Operator
        xjin2 = spg_bp(A,b,options); %Use basis pursuit to recover our solution
        djin2 = C'*xjin2; %Transform recovery back to time-space domain
        Djin2 = reshape(djin2,nr,ns); 
        recon_jt(i:i+k-1,j:j+k-1) = Djin2;
        Ddiffjin2 = D - Djin2; %Our Error
        SNRjin2 = mean(abs(Ddiffjin2(:)).^2);
        jit_err_2d = jit_err_2d + SNRjin2;  
    end
end

uni_err = uni_err/9
jit_err_1d = jit_err_1d/9
jit_err_2d = jit_err_2d/9


getd('C:\Users\arafat\Dropbox\Elastography');
env = envelope(D1);
env1 = envelope(recon_ur);
env_s = envelope(recon_jo);
env_ns = envelope(recon_jt);
env = medfilt2(env);
env1 = medfilt2(env1,[9,1]);
env_s= medfilt2(env_s,[9,1]);
env_ns = medfilt2(env_ns);

option = 'mu';
fac = 255;
y = LogComp(env,option,fac);
y = 1 + y./max(y(:));
y1 = LogComp(env1,option,fac);
y1 = 1 + y1./max(y1(:));
y2 = LogComp(env_s,option,fac);
y2 = 1 + y2./max(y2(:));
y3 = LogComp(env_ns,option,fac);
y3 = 1 + y3./max(y3(:));

MSE_un = sum(sum((y-y1).^2))/size(y1,1)/size(y1,2);
MSE_jo = sum(sum((y-y2).^2))/size(y2,1)/size(y2,2);
MSE_jt = sum(sum((y-y3).^2))/size(y3,1)/size(y3,2);

[ssim1, ssim_map] = ssim(y, y1, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM_uniform = ssim1
[ssim2, ssim_map] = ssim(y, y2, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM_jitter_1d = ssim2
[ssim3, ssim_map] = ssim(y, y3, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
MSSIM_jitter_2d = ssim3

figure; imagesc(y);colormap gray;axis off
figure; imagesc(y1);colormap gray;axis off
figure; imagesc(y2);colormap gray;axis off
figure; imagesc(y3);colormap gray;axis off