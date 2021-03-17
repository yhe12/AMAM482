%% clc;clear all
load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')
numFrames1_1 = size(vidFrames1_1,4); 
numFrames2_1 = size(vidFrames2_1,4); 
numFrames3_1 = size(vidFrames3_1,4); 
data1 = [];
for i = 1:numFrames1_1
    filter = zeros(480, 640);
    filter(150:250, 350:450) = 1;
    X = vidFrames1_1(:,:,:,i);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data1 = [data1; mean(X), mean(Y)];
end
 
data2 = [];
filter = zeros(480, 640);
filter(170:430, 400:500) = 1;
for j = 1:numFrames2_1
    X = vidFrames2_1(:,:,:,j);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data2 = [data2; mean(X), mean(Y)];
end
 
data3 = [];
filter = zeros(480, 640);
filter(200:350, 235:485) = 1;
for k = 1:numFrames3_1
    X = vidFrames3_1(:,:,:,k);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data3 = [data3; mean(X), mean(Y)];
end 
 
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);
data1 = data1(1:length(data3),:);
data2 = data2(1:length(data3),:);
alldata = [data1';,data2';data3']
[m,n] = size(alldata);
avg = mean(alldata,2);
all = alldata-repmat(avg,1,n);
[u,s,v]=svd(all/sqrt(n-1));
sig = diag(s);
figure(1)
plot(1:6, sig.^2/sum(sig.^2),'ko', 'Markersize', 10)
title("Test 1")
figure(2)
plot(1:214,all(2,:),1:214,all(4,:),1:214,all(6,:), 'Linewidth',2)
title("Test 1")
 
%% clc;clear all
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')
numFrames1_2 = size(vidFrames1_2,4); 
numFrames2_2 = size(vidFrames2_2,4); 
numFrames3_2 = size(vidFrames3_2,4); 
data1 = [];
for i = 1:numFrames1_2
    filter = zeros(480, 640);
    filter(150:250, 350:450) = 1;
    X = vidFrames1_2(:,:,:,i);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data1 = [data1; mean(X), mean(Y)];
end
 
data2 = [];
filter = zeros(480, 640);
filter(170:430, 400:500) = 1;
for j = 1:numFrames2_2
    X = vidFrames2_2(:,:,:,j);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data2 = [data2; mean(X), mean(Y)];
end
 
data3 = [];
filter = zeros(480, 640);
filter(200:350, 235:485) = 1;
for k = 1:numFrames3_2
    X = vidFrames3_2(:,:,:,k);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data3 = [data3; mean(X), mean(Y)];
end 
 
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);
data2 = data1(1:length(data1),:);
data3 = data2(1:length(data1),:);
alldata = [data1';,data2';data3']
[m,n] = size(alldata);
avg = mean(alldata,2);
all = alldata-repmat(avg,1,n);
[u,s,v]=svd(all/sqrt(n-1));
sig = diag(s);
figure(3)
plot(1:6, sig.^2/sum(sig.^2),'ko', 'Markersize', 10)
title("Test 2")
figure(4)
plot(1:297,all(2,:),1:297,all(4,:),1:297,all(6,:), 'Linewidth',2)
title("Test 2")
%% clc;clear all
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')
numFrames1_3 = size(vidFrames1_3,4); 
numFrames2_3 = size(vidFrames2_3,4); 
numFrames3_3 = size(vidFrames3_3,4); 
data1 = [];
for i = 1:numFrames1_3
    filter = zeros(480, 640);
    filter(150:250, 350:450) = 1;
    X = vidFrames1_3(:,:,:,i);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data1 = [data1; mean(X), mean(Y)];
end
 
data2 = [];
filter = zeros(480, 640);
filter(170:430, 400:500) = 1;
for j = 1:numFrames2_3
    X = vidFrames2_3(:,:,:,j);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data2 = [data2; mean(X), mean(Y)];
end
 
data3 = [];
filter = zeros(480, 640);
filter(200:350, 235:485) = 1;
for k = 1:numFrames3_3
    X = vidFrames3_3(:,:,:,k);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data3 = [data3; mean(X), mean(Y)];
end 
 
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);
data2 = data1(1:length(data1),:);
data3 = data2(1:length(data1),:);
alldata = [data1';,data2';data3']
[m,n] = size(alldata);
avg = mean(alldata,2);
all = alldata-repmat(avg,1,n);
[u,s,v]=svd(all/sqrt(n-1));
sig = diag(s);
figure(5)
plot(1:6, sig.^2/sum(sig.^2),'ko', 'Markersize', 10)
title("Test 3")
figure(6)
plot(1:237,all(2,:),1:237,all(4,:),1:237,all(6,:), 'Linewidth',2)
title("Test 3")
%% clc;clear all
load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')
numFrames1_4 = size(vidFrames1_4,4); 
numFrames2_4 = size(vidFrames2_4,4); 
numFrames3_4 = size(vidFrames3_4,4); 
data1 = [];
for i = 1:numFrames1_4
    filter = zeros(480, 640);
    filter(150:250, 350:450) = 1;
    X = vidFrames1_4(:,:,:,i);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data1 = [data1; mean(X), mean(Y)];
end
 
data2 = [];
filter = zeros(480, 640);
filter(170:430, 400:500) = 1;
for j = 1:numFrames2_4
    X = vidFrames2_4(:,:,:,j);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data2 = [data2; mean(X), mean(Y)];
end
 
data3 = [];
filter = zeros(480, 640);
filter(200:350, 235:485) = 1;
for k = 1:numFrames3_4
    X = vidFrames3_4(:,:,:,k);
    Xg = rgb2gray(X);
    X2 = double(X);
    Xg2 = double(Xg);
    Xf = Xg2.*filter;
    thres = Xf > 10;
    ind = find(thres);
    [Y, X] = ind2sub(size(thres), ind);
    data3 = [data3; mean(X), mean(Y)];
end 
 
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);
data2 = data1(1:length(data1),:);
data3 = data2(1:length(data1),:);
alldata = [data1';,data2';data3']
[m,n] = size(alldata);
avg = mean(alldata,2);
all = alldata-repmat(avg,1,n);
[u,s,v]=svd(all/sqrt(n-1));
sig = diag(s);
figure(7)
plot(1:6, sig.^2/sum(sig.^2),'ko', 'Markersize', 10)
title("Test 4")
figure(8)
plot(1:389,all(2,:),1:389,all(4,:),1:389,all(6,:), 'Linewidth',2)
title("Test 4")
