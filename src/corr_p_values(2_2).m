%% clear env,get and set current directory
clc;
clear;
currdir = pwd;
fprintf(currdir);
userpath(currdir); %set working directory to current dir of .m file

%% load data
%loaddata;
name = '/energydata_complete.csv';
filename = strcat(currdir,name)
data = importfile(filename)
disp('Hi0');
disp('Hi00');

%% Plot all variables with respect to time
figure;
stackedplot(data)
%% convert date to number 
data.date = datenum(data.date, 'yyyy-mm-dd HH:MM:SS');
ts = data.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data.date = ts;
disp('Hi');

%% set parameters
alpha = 0.05;
K = 29;
nonparametric = 1; 
rthresh = 0.2;
M = 100;    
maxwordlength = 30;
rng(1);
[m,n]=size(data);
disp('Hi2');
%% Select subset of K genes (K<=m)
iV = randperm(m);
data2 = data(iV(1:K),:);
tmpdata = table2array(data2)
disp('Hi3');
%% Compute the correlation matrix and the significance (p-values)
[ccM,p1M] = corrcoef(tmpdata); % compute the correlation coefficients and p-values of my matrix
p1M(1:K+1:K*K) = 0; % assign diagonal values to zero
if nonparametric
    % Randomized r values for the pair (x1,x2)
    p2M = zeros(K,K);
    ccsurT = NaN*ones(M,K,K);
    ccsurT(1,:,:) = ccM;
    for i=1:M
        zM = NaN*ones(n,K);
        for j=1:K
            zM(:,j) = tmpdata(randperm(n),j);
        end
        ccsurT(i+1,:,:) = corrcoef(zM);
    end
    for ik=1:K-1
        for jk=ik+1:K
            rxyV = squeeze(ccsurT(:,ik,jk));
            p2M(ik,jk) = resampledpvalue(rxyV,1);
            p2M(jk,ik) = p2M(ik,jk);
        end
    end
end
tit1txt = sprintf('R_{XY}');
C = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28', '29'}
h1 = plotnetworktitle(abs(ccM),[], C,tit1txt,1,0,maxwordlength);

adj1M = p1M < alpha;
tit2txt = sprintf('Parametric p(R_{XY}) < %1.2f',alpha);
h2 = plotnetworktitle(adj1M,[0 1],C,tit2txt,2,0,maxwordlength);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f), parametric R_{XY}',alpha);
h3 = plotnetworktitle(adjfdr1M,[0 1],C,tit3txt,3,0,maxwordlength);

adjthrM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY} > %1.2f',rthresh);
h4 = plotnetworktitle(adjthrM,[0 1],C,tit4txt,4,0,maxwordlength);

if nonparametric
    adj2M = p2M < alpha;
    tit5txt = sprintf('Randomization p(R_{XY}) < %1.2f',alpha);
    h5 = plotnetworktitle(adj2M,[0 1],C,tit5txt,5,0,maxwordlength);
    adjfdr2M = adjFDRmatrix(p2M,alpha,2);
    tit6txt = sprintf('FDR (a=%1.3f), randomization R_{XY}',alpha);
    h6 = plotnetworktitle(adjfdr2M,[0 1],C,tit6txt,6,0,maxwordlength);
end

disp('Hi4');
%% Explore autocorrelations and cross correlations in greek stocks
maxtau = 20;
maxtau2 = 10;
p = 1;
M = 100;
alpha = 0.05;
% zalpha = norminv(1-alpha/2);
zalpha = 1.96;

yM = load('../data/stocks2003.dat');
[n,m]=size(yM);
% rng(1);

% Read the names of the stocks
nameM = textread('../data/stock_names.dat','%s');