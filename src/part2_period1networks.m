% Lab 5: cross correlations in greek stocks and networks.
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% configuration
tau = 0;
alpha = 0.01;
K = 5;
rthresh = 0.05;
maxtau = 20;
p=5;
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
nameM = nameM(:,2:end);
%%nameM = nameM';
tmpdata = table2array(data1);
tmpdata = tmpdata(:,2:end); %remove 1st column (date)
[n,m]=size(tmpdata);
%yM = tmpdata';
yM=tmpdata;
K=m;


% rng(1);

%iV = randperm(m);
%yM = yM(:,iV(1:K));

% If NaN replace them with interpolated values for each time series
for i=1:K
    i1V = find(isnan(yM(:,i)));
    if ~isempty(i1V)
        iokV = setdiff([1:n]',i1V);
        yM(i1V,i) = interp1(iokV,yM(iokV,i),i1V,'spline');
    end
end

% Use log returns
xM=zeros(n,m-1);
for i=1:m
    y1=yM(:,i);
    
    xM(:,i)=fitAR(y1,p);
end

% For each pair compute the correlation matrix
ccM = NaN*ones(K,K);
p1M = zeros(K,K);

% replace NaN values to handle the problem of issymmetric() - no affect 
ccM(isnan(ccM)) = 0; 
isccMSymmetric = false;
% check if is symmetric
if issymmetric(ccM)
    disp("ccM is symmetric. Answer Y.");
    isccMSymmetric = true;
else
   disp("ccM is asymmetric. Answer N.");
end


adj1M = p1M < alpha;

% replace NaN values to handle the problem of issymmetric() - no affect 
p1M(isnan(p1M)) = 0; 
isp1MSymmetric = false;
% check if is symmetric
if issymmetric(p1M)
    disp("p1M is symmetric. Answer Y.");
    isp1MSymmetric = true;
else
   disp("p1M is asymmetric. Answer N.");
end

if tau==0
    % The correlation matrix is symmetric
    [ccM,p1M] = corrcoef(xM);
    p1M(1:K+1:K*K) = 0;
else
    % The correlation matrix is not symmetric
    for ik=1:K-1
        for jk=ik+1:K
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,ik),xM(1+tau:end,jk));
            ccM(ik,jk) = tmpM(1,2);
            p1M(ik,jk) = ptmpM(1,2);
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,jk),xM(1+tau:end,ik));
            ccM(jk,ik) = tmpM(1,2);
            p1M(jk,ik) = ptmpM(1,2);
        end
    end
end    
tit1txt = sprintf('R_{XY}(%d)',tau);
h1 = plotnetworktitle(ccM,[],nameM,tit1txt,1);



% replace NaN values to handle the problem of issymmetric() - no affect 
ccM(isnan(ccM)) = 0; 
isccMSymmetric = false;
% check if is symmetric
if issymmetric(ccM)
    disp("ccM is symmetric. Answer Y.");
    isccMSymmetric = true;
else
   disp("ccM is asymmetric. Answer N.");
end


adj1M = p1M < alpha;

% replace NaN values to handle the problem of issymmetric() - no affect 
p1M(isnan(p1M)) = 0; 
isp1MSymmetric = false;
% check if is symmetric
if issymmetric(p1M)
    disp("p1M is symmetric. Answer Y.");
    isp1MSymmetric = true;
else
   disp("p1M is asymmetric. Answer N.");
end

tit2txt = sprintf('Adjacency p(R_{XY}(%d)) < %1.2f',tau,alpha);
h2 = plotnetworktitle(adj1M,[0 1],nameM,tit2txt,2);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) R_{XY}(%d)',alpha,tau);
h3 = plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,3);

rthreshM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY}(%d) > %1.2f',tau,rthresh);
h4 = plotnetworktitle(rthreshM,[0 1],nameM,tit4txt,4);

