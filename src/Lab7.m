% Lab7: Granger causality index (GCI) in EEG or financial data
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
data(:,1) = []; %remove col with timestamps
data(:,end) = []; %remove ranodm var1
data(:,end) = []; %remove random var2
tmpdata = table2array(data);
%% parameter set
alpha = 0.05; % significance level
K = 5; % Number of variables (time series) to use from the set of variables read in.
P = 10; % The order of the VAR model used for the computation of the 
        % Granger causality index (GCI) 
GCIthresh = 0.1; 
taus = 1/100; % The sampling time
rng(1);
fignow = 0;

xM = tmpdata';
%xM = xM(1:1000,1:end-1); % The last channel is ECG
[n,m]=size(xM);
% Read the names of the channels
txtC = data.Properties.VariableNames;
txtC = txtC';
iV = randperm(m);
% iV = [1:m];
%xM = xM(:,iV(1:K));
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
% xM = log(xM(2:n,:))-log(xM(1:n-1,:));

%% Show the multivariate time series
plotmts(xM,1,0,K,taus,nameM,fignow+1);

%% For each pair of channels compute the Granger causality index (GCI) and 
% form the GCI-causality matrix
fprintf('Computes the GCI (p=%d) for all %d variables...\n',P,K);
[GCIM,pGCIM] = GCI(xM,P,1);

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
