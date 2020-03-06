% Conditional Granger causality index (CGCI) in EEG or financial data
% and networks
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

data.date = datenum(data.date, 'yyyy-mm-dd HH:MM:SS');
ts = data.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples
                 %(to start from zero secs)
data.date = ts;
disp('Hi');

%% set parameters
alpha = 0.05; % significance level
K = 5; % Number of variables (time series) to use from the set of variables read in.
P = 10; % The order of the VAR model used for the computation of the 
        % conditional Granger causality index (CGCI) 
        % var is vector auto regression=>https://en.wikipedia.org/wiki/Vector_autoregression
CGCIthresh = 0.1; 
taus = 600; % The sampling time
rng(1);
fignow = 5;

tmpdata = table2array(data);
tmpdata = tmpdata';
xM = tmpdata
[n,m]=size(xM);
% Read the names of the channels
iV = randperm(n);

txtC = data.Properties.VariableNames;
txtC = txtC';

xM = xM(:,iV(1:K));
nameM = txtC(iV(1:K),:);


% If NaN replace them with interpolated values for each time series
for i=1:K
    i1V = find(isnan(xM(:,i)));
    if ~isempty(i1V)
        iokV = setdiff([1:n]',i1V);
        xM(i1V,i) = interp1(iokV,xM(iokV,i),i1V,'spline');
    end
end

% Use log returns
xM = log(xM(2:n,:))-log(xM(1:n-1,:));

%% Show the multivariate time series
plotmts(xM,1,0,K,taus,nameM,fignow+1);

%% For each pair of channels compute the Granger causality index (GCI) and 
% form the GCI-causality matrix
fprintf('Computes the CGCI (p=%d) for all %d variables...\n',P,K);
[CGCIM,pCGCIM] = CGCI(xM,P,1);

%% Plot the CGCI-causality network
% The network of weighted connections given by CGCI_{X->Y}(P)
tit1txt = sprintf('CGCI_{X->Y}(%d)',P);
plotnetworktitle(CGCIM,[],nameM,tit1txt,fignow+2);

adj1M = pCGCIM < alpha;
tit2txt = sprintf('Adjacency p(CGCI_{X->Y}(%d)) < %1.2f',P,alpha);
plotnetworktitle(adj1M,[0 1],nameM,tit2txt,fignow+3);

adjfdr1M = adjFDRmatrix(pCGCIM,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) CGCI_{X->Y}(%d)',alpha,P);
plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,fignow+4);

CGCIthreshM = CGCIM > CGCIthresh;
tit4txt = sprintf('Adjacency CGCI_{X->Y}(%d) > %1.2f',P,CGCIthresh);
plotnetworktitle(CGCIthreshM,[0 1],nameM,tit4txt,fignow+5);
