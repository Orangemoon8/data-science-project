clear; close all; clc
% HW4
%% Case 1: Band Classification
%load music data
path = './song1/';
file = dir(path);
%iterate through each song
testmatrix = [];
for i = 3:length(file) %length(file) = 8
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 20:5:90
       clip = y(1,Fs*t:Fs*(t+5));
       clipsp = abs(spectrogram(clip));
       clipfft = reshape(clipsp,[1,8*16385]);
       testmatrix = [testmatrix; clipfft];
   end
end
testmatrix = testmatrix';
feature = 20;
[U,S,V,w,sortclassical,sortjazz,sortpop] = trainer(testmatrix,feature);
figure(1)
hold on
nc = size(testmatrix,2)/3; nj = size(testmatrix,2)/3; np = size(testmatrix,2)/3;
plot(sortclassical, zeros(nc), 'go')
plot(sortjazz, 0.5*ones(nj), 'ro')
plot(sortpop, ones(np), 'ko')
title('projection onto the w')
ylabel('artist/genre value')
xlabel('data point')
% From obervation of the graph, classical < threshold < jazz < pop
t1 = length(sortclassical);
t2 = 1;
while sortclassical(t1)>sortjazz(t2)
    t1 = t1-1;
    t2 = t2+1;
end
th1 = (sortclassical(t1)+sortjazz(t2))/2;
t3 = length(sortjazz);
t4 = 1;
while sortjazz(t3)>sortpop(t4)
    t3 = t3-1;
    t4 = t4+1;
end
th2 = (sortjazz(t3)+sortpop(t4))/2;
xline(th1,'--');
xline(th2, '.-');
print(gcf,'-dpng','figure 1.png');
% Test Classifier
%(1 classical, 1 jazz, 1 pop)
path = './test1/';
file = dir(path);
%iterate through each song
%(2 classical, 2 jazz, 2 pop) I may increase the number of songs downloaded
test = [];
for i = 3:length(file) %length(file) = 8
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 20:5:70
       tclip = y(1,Fs*t:Fs*(t+5));
       tclipsp = abs(spectrogram(tclip));
       tclipfft = reshape(tclipsp,[1,8*16385]);
       test = [test; tclipfft];
   end
end
test = test';

nc2 = size(test,2)/3; nj2 = size(test,2)/3; np2 = size(test,2)/3;
classicallabel = zeros(1,nc2);
jazzlabel = ones(1,nj2)*0.5;
poplabel = ones(1,np2);
hiddenlabel = [classicallabel, jazzlabel, poplabel];
TestMat = U'*test;  % PCA projection
pval = w'*TestMat;  % LDA projection

totaltest = nc2 + nj2 + np2;
record = [];
for i = 1: totaltest
   if pval(i) <= th1
       record = [record;0];
   elseif pval(i) <= th2
       record = [record;0.5];
   else 
       record = [record; 1];  
   end      
end
record = record';
errNum = sum(abs(record - hiddenlabel)~=0);
accuracy1 = 1 - errNum/totaltest

%% Case 2: The case for Seattle
%load music data
path = './song2/';
file = dir(path);
%iterate through each song
testmatrix = [];
for i = 3:length(file)
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 10:5:100
       clip = y(1,Fs*t:Fs*(t+5));
       clipsp = abs(spectrogram(clip));
       clipfft = reshape(clipsp,[1,8*16385]);
       testmatrix = [testmatrix; clipfft];
   end
end
testmatrix = testmatrix';
figure(2)
feature = 20;
[U,S,V,w,sortclassical,sortjazz,sortpop] = trainer(testmatrix,feature);
% (green)classical here is Derek's music
% (red)Jazz here is Lobo's music
% (black)Pop here is Thorn's music
hold on
nc = size(testmatrix,2)/3; nj = size(testmatrix,2)/3; np = size(testmatrix,2)/3;
plot(sortclassical, zeros(nc), 'go')
plot(sortjazz, 0.5*ones(nj), 'ro')
plot(sortpop, ones(np), 'ko')
title('projection onto the w')
ylabel('artist/genre value')
xlabel('data point')
% From obervation of the graph, classical < threshold < jazz < pop
t1 = length(sortjazz);
t2 = 1;
while sortjazz(t1)>sortclassical(t2)
    t1 = t1-1;
    t2 = t2+1;
end
th1 = (sortclassical(t1)+sortjazz(t2))/2;
t3 = length(sortclassical);
t4 = 1;
while sortclassical(t3)>sortpop(t4)
    t3 = t3-1;
    t4 = t4+1;
end
th2 = (sortclassical(t3)+sortpop(t4))/2;
xline(th1,'--');
xline(th2, '.-');
print(gcf,'-dpng','figure 2.png');
% Test Classifier
%(1 classical, 1 jazz, 1 pop)
path = './test2/';
file = dir(path);
%iterate through each song
test = [];
for i = 3:length(file) 
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 10:5:80
       tclip = y(1,Fs*t:Fs*(t+5));
       tclipsp = abs(spectrogram(tclip));
       tclipfft = reshape(tclipsp,[1,8*16385]);
       test = [test; tclipfft];
   end
end
test = test';

nc2 = size(test,2)/3; nj2 = size(test,2)/3; np2 = size(test,2)/3;
classicallabel = zeros(1,nc2);
jazzlabel = ones(1,nj2)*0.5;
poplabel = ones(1,np2);
hiddenlabel = [classicallabel, jazzlabel, poplabel];
TestMat = U'*test;  % PCA projection
pval = w'*TestMat;  % LDA projection

totaltest = nc2 + nj2 + np2;
record = [];
for i = 1: totaltest
   if pval(i) <= th1
       record = [record;0.5];
   elseif pval(i) <= th2
       record = [record;0];
   else 
       record = [record; 1];  
   end      
end
record = record';
errNum = sum(abs(record - hiddenlabel)~=0);
accuracy2 = 1 - errNum/totaltest
%% Case 3: Genre Classification
%load music data
path = './song3/';
file = dir(path);
%iterate through each song
testmatrix = [];
for i = 3:length(file)
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 30:5:70
       clip = y(1,Fs*t:Fs*(t+5));
       clipsp = abs(spectrogram(clip));
       clipfft = reshape(clipsp,[1,8*16385]);
       testmatrix = [testmatrix; clipfft];
   end
end
testmatrix = testmatrix';
figure(3)
feature = 20;
[U,S,V,w,sortclassical,sortjazz,sortpop] = trainer(testmatrix,feature);
% (green)classical here is Derek's music
% (red)Jazz here is Lobo's music
% (black)Pop here is Thorn's music
hold on
nc = size(testmatrix,2)/3; nj = size(testmatrix,2)/3; np = size(testmatrix,2)/3;
plot(sortclassical, zeros(nc), 'go')
plot(sortjazz, 0.5*ones(nj), 'ro')
plot(sortpop, ones(np), 'ko')
title('projection onto the w')
ylabel('artist/genre value')
xlabel('data point')
% From obervation of the graph, classical < threshold < jazz < pop
t1 = length(sortpop);
t2 = 1;
while sortpop(t1)>sortjazz(t2)
    t1 = t1-1;
    t2 = t2+1;
end
th1 = (sortjazz(t1)+sortpop(t2))/2;
t3 = length(sortjazz);
t4 = 1;
while sortjazz(t3)>sortclassical(t4)
    t3 = t3-1;
    t4 = t4+1;
end
th2 = (sortclassical(t3)+sortjazz(t4))/2;
xline(th1,'--');
xline(th2, '.-');
print(gcf,'-dpng','figure 3.png');
% Test Classifier
%(1 classical, 1 jazz, 1 pop)
path = './test3/';
file = dir(path);
%iterate through each song
%(2 classical, 2 jazz, 2 pop) I may increase the number of songs downloaded
test = [];
for i = 3:length(file) %length(file) = 8
   filename = file(i).name;
   [y,Fs] = audioread(strcat(path, '/', filename));
   % Average the signal , treating it as mono type
   y = (y(:,1) + y(:,2))./2;
   y = y';
   % Cut 5 second clips as a sample in each song, y
   y = resample(y,Fs/2,Fs);
   Fs = Fs/2;
   for t = 10:5:40
       tclip = y(1,Fs*t:Fs*(t+5));
       tclipsp = abs(spectrogram(tclip));
       tclipfft = reshape(tclipsp,[1,8*16385]);
       test = [test; tclipfft];
   end
end
test = test';

nc2 = size(test,2)/3; nj2 = size(test,2)/3; np2 = size(test,2)/3;
classicallabel = zeros(1,nc2);
jazzlabel = ones(1,nj2)*0.5;
poplabel = ones(1,np2);
hiddenlabel = [classicallabel, jazzlabel, poplabel];
TestMat = U'*test;  % PCA projection
pval = w'*TestMat;  % LDA projection

totaltest = nc2 + nj2 + np2;
record = [];
for i = 1: totaltest
   if pval(i) <= th1
       record = [record;1];
   elseif pval(i) <= th2
       record = [record;0.5];
   else 
       record = [record; 0];  
   end      
end
record = record';
errNum = sum(abs(record - hiddenlabel)~=0);
accuracy3 = 1 - errNum/totaltest
%%
function [U,S,V, w,sortclassical,sortjazz,sortpop] = trainer(testmatrix,feature)
    nc = size(testmatrix,2)/3; nj = size(testmatrix,2)/3; np = size(testmatrix,2)/3; 
    
    [U,S,V] = svd(testmatrix,'econ');
    
    artist = S*V'; % projection onto principal components
    U = U(:,1:feature);
    classical = artist(1:feature,1:nc);
    jazz = artist(1:feature,nc+1:nc+nj);
    pop = artist(1:feature,nc+nj+1:nc+nj+np);
    
    mc = mean(classical,2);
    mj = mean(jazz,2);
    mp = mean(pop,2);
    totalm = mean(artist(1:feature),2);
    
    Sw = 0; % within class variances
    for k=1:nc
        Sw = Sw + (classical(:,k)-mc).*(classical(:,k)-mc)';
    end
    for k=1:nj
        Sw = Sw + (jazz(:,k)-mj).*(jazz(:,k)-mj)';
    end
    for k=1:np
       Sw = Sw + (pop(:,k)-mp).*(pop(:,k)-mp); 
    end
    
    Sb = (mc-totalm)*(mc-totalm)'+(mj-totalm)*(mj-totalm)'+(mp-totalm)*(mp-totalm)'; % between class 
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    vclassical = w'*classical; % projection onto the threshold
    vjazz = w'*jazz;
    vpop = w'*pop;
    
    sortclassical = sort(vclassical);
    sortjazz = sort(vjazz);
    sortpop = sort(vpop);
end