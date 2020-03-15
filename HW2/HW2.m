% Homework 1
clear all; close all; clc

%% Part 1
load handel
v = y';
% plot((1:length(v))/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');
% p8 = audioplayer(v,Fs);
% playblocking(p8);

L=length(v)/Fs;     % music time length
v = v(1:end-1);   % periodic
n=length(v);
t2=linspace(0,L,n+1);
t=t2(1:n); 
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

% Gaussian Filter with fixed translation 0.1 and changing window width
figure(1)
a_vec = [60 20 2 0.2];
for jj = 1:length(a_vec)
    a = a_vec(jj);
    tslide=0:0.1:L;
    Sgt_spec = zeros(length(tslide),n);
    Sgt_spec = [];
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2); 
        Sg=g.*v; 
        Sgt=fft(Sg); 
        Sgt_spec(j,:) = fftshift(abs(Sgt)); 
    end
    
    subplot(2,2,jj)
    pcolor(tslide,ks,Sgt_spec.'), 
    shading interp 
    title(['a = ',num2str(a)],'Fontsize',12)
    xlabel('Time [sec]');
    ylabel('Frequency [\omega]');
    colormap(hot) 
end
print(gcf, '-dpng', 'figure 1.png');
% Gaussian filter with fixed width a = 70 and changing translation
figure(2)
translation_vec = [3 1 0.1 0.05];
for jj = 1:length(translation_vec)
    a = 70;
    translation = translation_vec(jj);
    tslide=0:translation:L;
    Sgt_spec = zeros(length(tslide),n);
    Sgt_spec = [];
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2); 
        Sg=g.*v; 
        Sgt=fft(Sg); 
        Sgt_spec(j,:) = fftshift(abs(Sgt)); 
    end
    
    subplot(2,2,jj)
    pcolor(tslide,ks,Sgt_spec.'), 
    shading interp 
    title(['translation = ',num2str(translation)],'Fontsize',12)
    xlabel('Time [sec]');
    ylabel('Frequency [\omega]');
    colormap(hot)
end
print(gcf, '-dpng', 'figure 2.png');

% Compare Gaussian, Mexican hat wavlet and Shannnon step-function window
figure(3)
a = 70;
tslide=0:0.1:L;
% Gaussian window with width a = 70
Sgt_spec = zeros(length(tslide),n);
Sgt_spec = [];
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    Sg=g.*v; 
    Sgt=fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end
subplot(3,1,1)
pcolor(tslide,ks,Sgt_spec.'), 
shading interp
title(['Gaussian window'])
xlabel('Time [sec]');
ylabel('Frequency [\omega]');
colormap(hot)

% Mexican hat wavelet with width a = 0.1
a = 0.1;
Sgt_spec = zeros(length(tslide),n);
Sgt_spec = [];
for j=1:length(tslide)
    g=2/(sqrt(3*a) * (pi)^(1/4)) * (1-((t-tslide(j))/a).^2).*exp(-(t-tslide(j)).^2 / (2*a^2)); 
    Sg=g.*v; 
    Sgt=fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end
subplot(3,1,2)
pcolor(tslide,ks,Sgt_spec.'), 
shading interp
title(['Mexican hat wavelet window'])
xlabel('Time [sec]');
ylabel('Frequency [\omega]');
colormap(hot)

% Shannon Step-function window with width a = 0.1
Sgt_spec = zeros(length(tslide),n);
Sgt_spec = [];
for j=1:length(tslide)
    g= (abs(t - tslide(j)) < a); 
    Sg=g.*v; 
    Sgt=fft(Sg); 
    Sgt_spec(j,:) = fftshift(abs(Sgt)); 
end
subplot(3,1,3)
pcolor(tslide,ks,Sgt_spec.'), 
shading interp
title(['Shannon Step-function window'])
xlabel('Time [sec]');
ylabel('Frequency [\omega]');
colormap(hot)
print(gcf, '-dpng', 'figure 3.png');
%% Part 2
close all; clear all;
% Piano
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (piano)');
% p8 = audioplayer(y,Fs); playblocking(p8);
y = y'/2;
Lp=16;     % music time length
n=length(y);
t=(1:length(y))/Fs; 
k=(2*pi/Lp)*[0:n/2-1 -n/2:-1]; 
ks_p=fftshift(k);
tslide_p = 0:0.1:Lp;
pgt_spec = [];
score_p = [];
a = 50;
for j=1:length(tslide_p)
    g=exp(-a*(t-tslide_p(j)).^2); 
    Sg=g.*y; 
    Sgt=fft(Sg); 
    pgt_spec = [pgt_spec; abs(fftshift(Sgt))/max(abs(Sgt))]; 
    [M, I] = max(abs(Sgt));
    score_p = [score_p; abs(k(I))];
end
figure(4)
% subplot(2,1,1)
% plot(tslide_p,score_p/(2*pi));
% xlabel('Time [sec]');
% ylabel('Frequency(Hz)');
% title('Center Frequency of Piano');
% subplot(2,1,2)
pcolor(tslide_p,(ks_p/(2*pi)),pgt_spec.')
shading interp
title(['Piano'])
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
ylim([100 1000])
colormap(hot)
drawnow
hold on
print(gcf, '-dpng', 'figure 4.png');
%% Recorder
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (recorder)');
% p8 = audioplayer(y,Fs); playblocking(p8);
y = y'/2;
Lr=14;     % music time length
n=length(y);
t=(1:length(y))/Fs; 
k=(2*pi/Lr)*[0:n/2-1 -n/2:-1]; 
ks_r=fftshift(k);
tslide_r = 0:0.1:Lr;
rgt_spec = [];
score_r = [];
a = 50;
for j=1:length(tslide_r)
    g=exp(-a*(t-tslide_r(j)).^2); 
    Sg=g.*y; 
    Sgt=fft(Sg); 
    rgt_spec = [rgt_spec; abs(fftshift(Sgt))/max(abs(Sgt))]; 
    [M, I] = max(abs((Sgt)));
    score_r = [score_r; abs(k(I))];
end
figure(5)
% subplot(2,1,1)
% plot(tslide_r,score_r/(2*pi));
% xlabel('Time [sec]');
% ylabel('Frequency(Hz)');
% title('Center Frequency of Recorder');
% subplot(2,1,2)
pcolor(tslide_r,ks_r/(2*pi),rgt_spec.')
shading interp
title(['Recorder'])
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
colormap(hot)
ylim([100 1200])
drawnow
hold on
print(gcf, '-dpng', 'figure 5.png');
