clc;
clear all;
close all;
[a, fs] = audioread('audioclip-1586598895-3757.mp4');
sound(a,fs);

plot(a);
title("Input Audio Signal");
figure;
x = a';
% x = a(:,1)';
%% converting to mono channel
% x = a(:,1)'+a(:,2)';
% peak = max(abs(x));
% x = x/peak;
% peakL = max(abs(a(:,1)));
% peakR = max(abs(a(:,2)));
% max_ch = max([peakL peakR]);
% x = x*max_ch;
%%----------------
Nstart=length(x);
if(rem(length(x),2)~=0), x=[x 0]; end
N=length(x);
nl=2:N-1;
xx=zeros(1,N);
%%---hilbert transform
N=length(x); Nh=ceil(N/2);
k=0:N-1;

H=-1i.*sign(Nh-k).*sign(k);
x_hilb=ifft( fft(x).*H );
h = x_hilb;

%%--------hilbert completed

%%non-negative energy operator, envelope of derivative
envl(nl)=(x(nl+1).^2 + x(nl-1).^2 + h(nl+1).^2 + h(nl-1).^2)./4 -(x(nl+1).*x(nl-1) + h(nl+1).*h(nl-1))./2;
envl_deri=[0 0 envl(3:end-2) 0 0];


envl_deri=envl_deri(1:Nstart-1);

%plot:
envl_deri=(envl_deri-min(envl_deri))/(max(envl_deri)-min(envl_deri));
subplot(211); plot(x); ylabel('amplitude');
subplot(212); plot(envl_deri,'-'); %ylim([0 0.5])
ylabel('Energy using non--ve freq EO');
figure; 
sp_enh = envl_deri;
%% Enhancing the signal
samples = 50*fs/1000;
smooth_sig = smoothdata(envl_deri,'movmean',ceil(samples));
plot(smooth_sig);
title("Smoothened signal");
figure;
smooth_sig = rescale(smooth_sig,0,0.5);
%% fod %%
sp_enh = diff(smooth_sig);
sp_enh=(sp_enh-min(sp_enh))/(max(sp_enh)-min(sp_enh));
sp_enh = 1.*sp_enh;
plot(sp_enh);
% title("Enhanced Energy Signal Contour wo smooth");
% figure;
samples = 50*fs/1000;
sp_enh = smooth(sp_enh,samples/N);
sp_enh=(sp_enh-min(sp_enh))/(max(sp_enh)-min(sp_enh));
sp_enh = sp_enh - mean(sp_enh(:));
pts  = find(sp_enh(1:end-1)>0 & sp_enh(2:end) < 0);
plot(sp_enh);
title("Enhanced Energy Signal Contour");
figure;
%% Finding Evidence plot
g = gausswin(1600,2); % 5 = window width / 2
g = diff(g);
y = conv(g,smooth_sig);
samples = 50*fs/1000;
y = smooth(y,samples);

subplot(411);plot(sp_enh);title("enhanced signal");
subplot(412);plot(smooth_sig);title("smooth signal");
subplot(413);plot(y);
title("Spectral +ve peaks and -ve peaks representing VOP and VEP resp ");
hold on
[~,idp] = findpeaks(y);
[~,idn] = findpeaks(-y);
mean_posi = mean(y(idp)>0);
mean_negi = mean(y(idn)<0);

j=1;

idp_upd=[];
idn_upd=[];
for i = 1:length(idp)
    if(y(idp(i))>0.003)
        idp_upd(j) =idp(i);
        j=j+1;
    end
    
end
j=1;
for i = 1:length(idn)
    if(y(idn(i))<-0.001)
        idn_upd(j) =idn(i);
        j=j+1;
    end
    
end
l = size(y);
sq_wave = zeros(l(1,1),1);
if(idp_upd(1)<idn_upd(1))
%        fprintf('---%i    %i\n',idp_upd(1),idn_upd(1));
        sq_wave(idp_upd(1):idn_upd(1)) = 0.5;
end    
i=2;
j=2;
while(i~=length(idp_upd)+1)
    if(idp_upd(i)<idn_upd(j) && (idp_upd(i)>idn_upd(j-1)))
%        fprintf('---%i    %i\n',idp_upd(i),idn_upd(j));
        sq_wave(idp_upd(i):idn_upd(j)) = 0.5;
        j=j+1;
        i=i+1;
    elseif (idp_upd(i)>idn_upd(j) && (idp_upd(i)>idn_upd(j-1)))
        j=j+1;
    elseif(idp_upd(i)<idn_upd(j-1))
        i=i+1;
    
    end
end
plot(sq_wave);ylim([-0.3 0.8]);
%title('Vowel Region');
hold off
subplot(414);plot(x);ylim([-0.3 0.8]);
hold on
plot(sq_wave);
title('Detected Vowel Region');
hold off
figure;

plot(y);ylim([-.3 0.8]);xlim([0 Nstart]);
title("Vowel Region Detection");
hold on
plot(sq_wave);
hold off
figure;

plot(x);ylim([-.3 0.8]);xlim([0 Nstart]);
title("Vowel Region Detection");
hold on
plot(sq_wave);
hold off
