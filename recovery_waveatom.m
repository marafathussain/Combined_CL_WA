% Load Data
close all;
clear all;
clc;

getd = @(p)path(p,path);
getd('C:\Users\arafat\Dropbox\WA\spot-1.2.wa');
getd('C:\Users\arafat\Dropbox\WA\spgl1-1.8.wa');
getd('C:\Users\arafat\Dropbox\WA\WaveAtom-1.1.1');
getd('C:\Users\arafat\Dropbox\CS US');

% load 'C:\Users\arafat\Dropbox\CS US\cyst_d40x40_c5s40_nrf.mat'
% D1 = 256.*(rf./max(rf(:)));
% D1 = downsample(D1,2);

% load 'C:\Users\arafat\Dropbox\Elastography\AM2D\rf052.mat'
% rf = RfDataFilt(1:2048,233:360);clear RfDataFilt;
% rf = downsample(rf,2);
% D1 = 256.*(rf./max(rf(:)));

% load 'C:\Users\arafat\Dropbox\Elastography\AM2D\rf021.mat'
% rf = RfDataFilt(1:2048,76:203);clear RfDataFilt;
% rf = downsample(rf,2);
% D1 = 256.*(rf./max(rf(:)));

% load 'C:\Users\arafat\Dropbox\Elastography\AM2D\rf011.mat'
% rf = RfDataFilt(1:2304,251:378);clear RfDataFilt;
% rf = downsample(rf,2);
% D1 = 256.*(rf./max(rf(:)));

load 'C:\Users\arafat\Dropbox\Elastography\AM2D\rf043.mat'
rf = RfDataFilt(1:2048,73:200);clear RfDataFilt;
rf = downsample(rf,2);
D1 = 256.*(rf./max(rf(:)));

recon_c = zeros(size(D1));
recon_w = recon_c;

[xx,yy] = size(D1);

% Define undersampling percentage
ss = [10 5 3 2];
ts = 1;
% perc = 1/ss;
% perct = 1/ts;
k = 128;

for mm = 2:2
    for nn = 1:1
        for i = 1:k:xx-k+1
            for j = 1:k:yy-k+1

                urm = zeros(k);

                D = D1(i:i+k-1,j:j+k-1);
                [nr,ns] = size(D);

                % Uniform random sampling along scan-line
                for g = 1:k
                    p = randperm(nr);
                    r = sort(p(1:floor(nr/(ss(mm)*ts))));
                    urm(r,g) = 1;
                end
                Runf = opMask(ns*nr,find(urm(:))); % Restriction Operator


                % Instruct spgl1 recovery with Curvelet
                CC = opCurvelet(nr,ns); %Sparse Transform


                % Instruct spgl1 recovery with Waveatom
                CW = opWaveatom(nr,ns); %Sparse Transfor


                % Instruct spgl1 recovery with Curvelet+Waveatom
                C = opCW(nr,ns); %Sparse Transfor


                options = spgSetParms('optTol', 5e-4, 'bpTol', 5e-4, ...
                                      'iterations', 100, 'verbosity', 0);

                % Curvelet Recovery
                RD = Runf*D(:);
                b = RD(:); % Observed subsampled data
                A = Runf*CC'; % Measurement Operator
                xunf = spg_bp(A,b,options); %Use basis pursuit to recover our solution
                dunf = CC'*xunf; %Transform recovery back to time-space domain
                Dunf = reshape(dunf,nr,ns); 
                recon_c(i:i+k-1,j:j+k-1) = Dunf;       

                % Waveatom Recovery
                RD = Runf*D(:);
                b = RD(:); % Observed subsampled data
                A = Runf*CW'; % Measurement Operator
                xunf = spg_bp(A,b,options); %Use basis pursuit to recover our solution
                dunf = CW'*xunf; %Transform recovery back to time-space domain
                Dunf = reshape(dunf,nr,ns); 
                recon_w(i:i+k-1,j:j+k-1) = Dunf;  

                % Combined Recovery
                RD = Runf*D(:);
                b = RD(:); % Observed subsampled data
                A = Runf*C'; % Measurement Operator
                xunf = spg_bp(A,b,options); %Use basis pursuit to recover our solution
                dunf = C'*xunf; %Transform recovery back to time-space domain
                Dunf = reshape(dunf,nr,ns); 
                recon(i:i+k-1,j:j+k-1) = Dunf; 
            end
        end
        
        if mm == 1 && nn == 1
            err_c = zeros(3,4);err_w = zeros(3,4);err = zeros(3,4);
            MSE_c = zeros(3,4);MSE_w = zeros(3,4);MSE = zeros(3,4);
            MSSIM_c = zeros(3,4);MSSIM_w = zeros(3,4);MSSIM = zeros(3,4);
        end

        %% Evaluation
        err_c(nn,mm) = mean(abs(D1(:) - recon_c(:)).^2);
        err_w(nn,mm) = mean(abs(D1(:) - recon_w(:)).^2);
        err(nn,mm)   = mean(abs(D1(:) - recon(:)).^2);

        getd('C:\Users\arafat\Dropbox\Elastography');

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


        MSE_c(nn,mm) = sum(sum((y-y1).^2))/size(y1,1)/size(y1,2);
        MSE_w(nn,mm) = sum(sum((y-y2).^2))/size(y2,1)/size(y2,2);
        MSE(nn,mm) = sum(sum((y-y3).^2))/size(y3,1)/size(y3,2);


        [ssim1, ssim_map] = ssim(y, y1, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
        MSSIM_c(nn,mm) = ssim1;
        [ssim2, ssim_map] = ssim(y, y2, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
        MSSIM_w(nn,mm) = ssim2;
        [ssim3, ssim_map] = ssim(y, y3, [0.01 0.03], fspecial('gaussian', 11, 1.5), 2);
        MSSIM(nn,mm) = ssim3;
    end
end

figure; imagesc(y);colormap gray;axis off;title('Actual')
figure; imagesc(y1);colormap gray;axis off;title('Curvelet')
figure; imagesc(y2);colormap gray;axis off;title('Waveatom')
figure; imagesc(y3);colormap gray;axis off;title('Combined')