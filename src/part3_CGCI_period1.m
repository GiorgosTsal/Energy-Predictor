%% Conditional Granger causality index (CGCI)
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% load my dataset and convert date to number 
name = '/energydata_complete.csv';
filename = strcat(currdir,name)
data = importfile(filename)
data1=data(1:6*24*63, :);
data2=data(6*24*63+1:6*24*63*2,:);

data1.date = datenum(data1.date, 'yyyy-mm-dd HH:MM:SS');
ts = data1.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data1.date = ts;

tmpdata = table2array(data1)
tmpdata = tmpdata(:,2:end);
%% set parameters
alpha = 0.01; % significance level
K = 26; % Number of variables (time series) to use from the set of variables read in.
P = 10; % The order of the VAR model used for the computation of the 
        % conditional Granger causality index (CGCI) 
        % var is vector auto regression=>https://en.wikipedia.org/wiki/Vector_autoregression
CGCIthresh = 0.01; 
taus = 600; % The sampling time
rng(1);
fignow = 1;

xM = tmpdata;
[n,m]=size(xM);

%assign variable names
nameM = data.Properties.VariableNames;
nameM = nameM(:,2:end);
%%
for i=1:m
    i1V = find(isnan(xM(:,i)));
    if ~isempty(i1V)
        iokV = setdiff([1:n]',i1V);
        xM(i1V,i) = interp1(iokV,xM(iokV,i),i1V,'spline');
    end
end
fprintf('Computes the CGCI (p=%d) for all %d variables...\n',P,K);
[CGCIM,pCGCIM] = CGCI(xM,P,1);
%% Plot the CGCI-causality network
% The network of weighted connections given by CGCI_{X->Y}(P)
tit1txt = sprintf('CGCI_{X->Y}(%d)',P);
%% check for symmetry
% replace NaN values to handle the problem of issymmetric() - no affect 
CGCIM(isnan(CGCIM)) = 0; 
isGCIMSymmetric = false;
% check if is symmetric
if issymmetric(CGCIM)
    disp("CGCIM is symmetric. Answer Y.");
    isCGCIMSymmetric = true;
else
   disp("CGCIM is asymmetric. Answer N.");
end

% replace NaN values to handle the problem of issymmetric() - no affect 
pCGCIM(isnan(pCGCIM)) = 0; 
ispCGCIMSymmetric = false;
% check if is symmetric
if issymmetric(pCGCIM)
    disp("pCGCIM is symmetric. Answer Y.");
    ispCGCIMSymmetric = true;
else
   disp("pCGCIM is asymmetric. Answer N.");
end
%% plot networks
plotnetworktitle(CGCIM,[],nameM,tit1txt,fignow);

adj1M = pCGCIM < alpha;
tit2txt = sprintf('Adjacency p(CGCI_{X->Y}(%d)) < %1.2f',P,alpha);
plotnetworktitle(adj1M,[0 1],nameM,tit2txt,fignow+1);

adjfdr1M = adjFDRmatrix(pCGCIM,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) CGCI_{X->Y}(%d)',alpha,P);
plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,fignow+2);

CGCIthreshM = CGCIM > CGCIthresh;
tit4txt = sprintf('Adjacency CGCI_{X->Y}(%d) > %1.2f',P,CGCIthresh);
plotnetworktitle(CGCIthreshM,[0 1],nameM,tit4txt,fignow+3);

