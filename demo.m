clear all
clc
close all
addpath('Database')
load('orl_56_46.mat')
 p =56;
 q =46;
ClassNum = length(unique(sample_label));
 sample_data= sample_data./255;
 train_num=3;
experiments=10;
reg_rate=zeros(1, experiments);
for ii = 1:experiments
    train_data = []; test_data = [];
    train_label = []; test_label = [];
    for i = 1 : ClassNum
        index = find(sample_label == i);
        randindex = index(randperm((length(index))));
        train_data = [train_data sample_data(:,randindex(1 : train_num))];
        train_label = [train_label  sample_label(randindex(1 : train_num))];
        test_data = [test_data sample_data(:,randindex(train_num + 1 : end))];
        test_label = [test_label  sample_label(randindex(train_num + 1 : end))];
    end  
ttdat=test_data;
trdat=train_data;
trls = train_label;
ttls= test_label;
imgsize = [p q];
n = size(ttdat, 2);
Proj = pinv(trdat'*trdat+0.01*eye(size(trdat,2)))*trdat';
alpha = 0.0001; 
delta=2;
ID=[];
for i = 1:n
    if mod(i,200)==0
        fprintf('%d / %d \n',i,n);
    end
    y = ttdat(:,i);
    Xs = Proj*y;
    Es = y-trdat*Xs;
    [w] = Weight(y, trdat, Xs, trls, delta);
   [X] = LDMR(y, trdat, w, trls, Xs, Es, alpha,  imgsize);
   [label] = classifier(trdat, X, trls, imgsize);  
    ID= [ID label];
end
 reg_rate(ii) = mean(ttls(:)==ID(:))
end
mean(reg_rate)


