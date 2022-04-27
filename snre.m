%%

load 'C:\Users\arafat\Dropbox\WA\EPstr_curvelet.mat';
strc = EPstr;clear EPstr;


a = strc(7:18,12:37); 
a1 = reshape(a,1,size(a,1)*size(a,2));
b = strc(7:18,90:115);
b1 = reshape(b,1,size(b,1)*size(b,2));
str3 = [a1 b1];
SNRe_c = mean(str3)/std(str3)

%%
load 'C:\Users\arafat\Dropbox\WA\EPstr_waveatom.mat';
strw = EPstr;clear EPstr;


a = strw(7:18,12:37); 
a1 = reshape(a,1,size(a,1)*size(a,2));
b = strw(7:18,90:115);
b1 = reshape(b,1,size(b,1)*size(b,2));
str3 = [a1 b1];
SNRe_w = mean(str3)/std(str3)


%%
load 'C:\Users\arafat\Dropbox\WA\EPstr_combined.mat';
str = EPstr;clear EPstr;


a = str(7:18,12:37); 
a1 = reshape(a,1,size(a,1)*size(a,2));
b = str(7:18,90:115);
b1 = reshape(b,1,size(b,1)*size(b,2));
str3 = [a1 b1];
SNRe = mean(str3)/std(str3)