%% Granger causality index (GCI)
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% load my dataset and convert date to number 
name = '/energydata_complete.csv';
filename = strcat(currdir,name);
data = importfile(filename);
data1=data(1:6*24*63, :);
data2=data(6*24*63+1:6*24*63*2,:);

data1.date = datenum(data1.date, 'yyyy-mm-dd HH:MM:SS');
ts = data1.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data1.date = ts;

%assign variable names
nameM = data1.Properties.VariableNames;
%%nameM = nameM';
tmpdata = table2array(data1)
tmpdata = tmpdata(:,2:end);
%% set parameters
alpha = 0.01; % significance level
K = 26; % Number of variables (time series) to use from the set of variables read in.
P = 10; % The order of the VAR model used for the computation of the 
        % conditional Granger causality index (CGCI) 
        % var is vector auto regression=>https://en.wikipedia.org/wiki/Vector_autoregression
GCIthresh = 0.01; 
taus = 600; % The sampling time
rng(1);
fignow = 5;

xM = tmpdata;
[n,m]=size(xM);
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

fprintf('Computes the GCI (p=%d) for all %d variables...\n',P,K);
[GCIM,pGCIM] = GCI(xM,P,1);
%% check for symmetry
% replace NaN values to handle the problem of issymmetric() - no affect 
GCIM(isnan(GCIM)) = 0; 
isGCIMSymmetric = false;
% check if is symmetric
if issymmetric(GCIM)
    disp("GCIM is symmetric. Answer Y.");
    isGCIMSymmetric = true;
else
   disp("GCIM is asymmetric. Answer N.");
end

% replace NaN values to handle the problem of issymmetric() - no affect 
pGCIM(isnan(pGCIM)) = 0; 
ispGCIMSymmetric = false;
% check if is symmetric
if issymmetric(pGCIM)
    disp("pGCIM is symmetric. Answer Y.");
    ispGCIMSymmetric = true;
else
   disp("pGCIM is asymmetric. Answer N.");
end
%% Plot the GCI-causality network
% The network of weighted connections given by GCI_{X->Y}(P)
tit1txt = sprintf('GCI_{X->Y}(%d)',P);
plotnetworktitle(GCIM,[],nameM,tit1txt,fignow+2);

adj1M = pGCIM < alpha;
tit2txt = sprintf('Adjacency p(GCI_{X->Y}(%d)) < %1.2f',P,alpha);
plotnetworktitle(adj1M,[0 1],nameM,tit2txt,fignow+3);

adjfdr1M = adjFDRmatrix(pGCIM,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) GCI_{X->Y}(%d)',alpha,P);
plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,fignow+4);

GCIthreshM = GCIM > GCIthresh;
tit4txt = sprintf('Adjacency GCI_{X->Y}(%d) > %1.2f',P,GCIthresh);
plotnetworktitle(GCIthreshM,[0 1],nameM,tit4txt,fignow+5);
