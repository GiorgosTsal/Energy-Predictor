% This example shows dimension reduction in low and high-dimensional
% regression (discards unnecessary predictors). 
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
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data.date = ts;
disp('Hi');
% xTrainM = table2array(data);
% xTrainM(:,1) = []
% xTrainM(:,end) = []
% xTrainM(:,end) = []
% [n, p] = size(xTrainM)
tmpdata = table2array(data);
tmpdata(:,1) = [];
tmpdata(:,end) = [];
tmpdata(:,end) = [];

%% Split into train and test set
% Cross validation (train: 70%, test: 30%)
cv = cvpartition(size(tmpdata,1),'HoldOut',0.3);
idx = cv.test;

% Separate to training and test data
dataTrain = tmpdata(~idx,:);
dataTest  = tmpdata(idx,:);
disp('hi');

% set target 
yTrainM=dataTrain(:,1); % target energy from appliances
dataTrain(:,1)=[];
xTrainM=dataTrain;
[n, p] = size(xTrainM);
%% Generate response data Y = X * beta + eps , where beta has just a
% number dtrue of nonzero components, and the noise eps is normal.
dtrue = 1;
%iV = randperm(p);
iV = (1:p)';
betaV = zeros(p,1);
betaV(iV(1:dtrue)) = (2*unidrnd(2,dtrue,1)-3).*(round(4*rand(dtrue,1))+1);
yV = yTrainM;
d = dtrue; % The dimension reduction.

TSS = sum((yV-mean(yV)).^2);
mxV = mean(xTrainM);
xcM = xTrainM - repmat(mxV,n,1); % centered data matrix
my = mean(yV);
ycV = yV - my;

[uM,sigmaM,vM] = svd(xcM,'econ');
r = size(sigmaM,1);

%% OLS  
bOLSV = vM * inv(sigmaM) * uM'* ycV;
% yfitOLSV = xcM * bOLSV + my; 
bOLSV = [my - mxV*bOLSV; bOLSV];
yfitOLSV = [ones(n,1) xTrainM] * bOLSV; 
resOLSV = yV-yfitOLSV; 
RSSOLS = sum(resOLSV.^2);
rsquaredOLS = 1 - RSSOLS/TSS;
figure(1)
clf
plot(yV,yfitOLSV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('OLS R^2=%1.4f',rsquaredOLS))
figure(2)
clf
plot(yV,resOLSV/std(resOLSV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('OLS')

%% Stepwise fit
[bstepV,sdbV,pvalV,inmodel,stats]=stepwisefit(xTrainM,yV);
bstep0 = stats.intercept;
bstepV = [bstep0; bstepV].*[1 inmodel]';
yfitstepV = [ones(n,1) xTrainM] * bstepV; 
resstepV = yV-yfitstepV; 
RSSstep = sum(resstepV.^2);
rsquaredstep = 1 - RSSstep/TSS;
figure(3)
clf
plot(yV,yfitstepV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('stepwise regression R^2=%1.4f',rsquaredstep))
figure(4)
clf
plot(yV,resstepV/std(resstepV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('stepwise regression')

%% PCR
lambdaV = zeros(r,1);
lambdaV(1:d) = 1;
bPCRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
bPCRV = [my - mxV*bPCRV; bPCRV];
yfitPCRV = [ones(n,1) xTrainM] * bPCRV; 
resPCRV = yfitPCRV - yV;     % Calculate residuals
RSSPCR = sum(resPCRV.^2);
rsquaredPCR = 1 - RSSPCR/TSS;
figure(5)
clf
plot(yV,yfitPCRV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('PCR R^2=%1.4f',rsquaredPCR))
figure(6)
clf
plot(yV,resPCRV/std(resPCRV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('PCR')

%% PLS
[Xloadings,Yloadings,Xscores,Yscores,bPLSV] = plsregress(xTrainM,yV,d);
yfitPLSV = [ones(n,1) xTrainM]*bPLSV;
resPLSV = yfitPLSV - yV;     % Calculate residuals
RSSPLS = sum(resPLSV.^2);
rsquaredPLS = 1 - RSSPLS/TSS;
figure(7)
clf
plot(yV,yfitPLSV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('PLS R^2=%1.4f',rsquaredPLS))
figure(8)
clf
plot(yV,resPLSV/std(resPLSV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('PLS')

%% Ridge regression
% [u2M,sigma2M,v2M] = svd(xcM);
% mu = (1/(n-p)) * sum((u2M(:,p+1:n)'*ycV).^2);
mu = RSSOLS/(n-p);
sigmaV = diag(sigmaM);
lambdaV = sigmaV.^2 ./ (sigmaV.^2 + mu);
bRRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
bRRV = [my - mxV*bRRV; bRRV];
% bRRV = ridge(yV,xTrainM,mu,0);
yfitRRV = [ones(n,1) xTrainM] * bRRV; 
resRRV = yfitRRV - yV;     % Calculate residuals
RSSRR = sum(resRRV.^2);
rsquaredRR = 1 - RSSRR/TSS;
figure(9)
clf
plot(yV,yfitRRV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('RR R^2=%1.4f',rsquaredRR))
figure(10)
clf
plot(yV,resRRV/std(resRRV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('RR')

%% LASSO 
[bM,fitinfo] = lasso(xcM,ycV);
lassoPlot(bM,fitinfo,'PlotType','Lambda','XScale','log');
% lambda = input('Give lambda >');
lambda = 0.5;
[lmin, ilmin] = min(abs(fitinfo.Lambda - lambda));
bLASSOV = bM(:,ilmin);
bLASSOV = [my - mxV*bLASSOV; bLASSOV];
yfitLASSOV = [ones(n,1) xTrainM] * bLASSOV; 
resLASSOV = yfitLASSOV - yV;     % Calculate residuals
RSSLASSO = sum(resLASSOV.^2);
rsquaredLASSO = 1 - RSSLASSO/TSS;
figure(11)
clf
plot(yV,yfitLASSOV,'.')
hold on
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('LASSO R^2=%1.4f',rsquaredLASSO))
figure(12)
clf
plot(yV,resLASSOV/std(resLASSOV),'.','Markersize',10)
hold on
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('LASSO')

fprintf('\t beta \t OLS \t step \t PCR \t PLS \t RR \t LASSO \n');
fprintf('\t Regression coefficient vectors \n');
disp([[0;betaV] bOLSV bstepV bPCRV bPLSV bRRV bLASSOV])
fprintf('\t R^2 values \n');
disp([NaN rsquaredstep rsquaredOLS rsquaredPCR rsquaredPLS rsquaredRR rsquaredLASSO])
