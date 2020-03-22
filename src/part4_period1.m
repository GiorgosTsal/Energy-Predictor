% This example shows dimension reduction in low and high-dimensional
% regression (discards unnecessary predictors). After estimating the
% regression model it uses it to make predictions in a test data matrix X
%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% load my dataset and convert date to number

d=5;
targetIdx=8;
name = '/energydata_complete.csv';
filename = strcat(currdir,name)


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

disp('Hi');
p = 10;
xM1 = table2array(data1);
xM1(:,1) = []
[n,m] = size(xM1);

trainIdx=round(0.7*n);
trainData=xM1(1:trainIdx,:);

testData=xM1(trainIdx+1:end,:);
yTrain=trainData(:,targetIdx);
yV=yTrain;

yTest=testData(:,targetIdx);
trainData(:,targetIdx)=[];
testData(:,targetIdx)=[];

[nTrain,m] = size(trainData);
[nTest,m] = size(testData);

xM=trainData;
ntest=nTest;
xtestM=testData;
ytestV=yTest;



TSStrain = sum((yTrain-mean(yTrain)).^2);
meanTrainData = mean(trainData);
mxV=meanTrainData;
centerTrainData = trainData - repmat(meanTrainData,nTrain,1); % centered data matrix
mYTrain = mean(yTrain);
centerYTrain = yTrain - mYTrain;

TSS=TSStrain;

TSStest = sum((yTest-mean(yTest)).^2);
meanTestData = mean(testData);
centerTestData = testData - repmat(meanTestData,nTest,1); % centered data matrix
my = mean(yTest);
centerYTest = yTest - my;

xcM=centerTrainData;
ycV=centerYTrain;

xctestM=centerTestData;
yctestV=centerYTest;

n=nTrain;

[uM,sigmaM,vM] = svd(xcM,'econ');
r = size(sigmaM,1);

%% OLS  
bOLSV = vM * inv(sigmaM) * uM'* ycV;
% yfitOLSV = xcM * bOLSV + my; 
bOLSV = [my - mxV*bOLSV; bOLSV];
yfitOLSV = [ones(n,1) xM] * bOLSV; 
resOLSV = yV-yfitOLSV; 
RSSOLS = sum(resOLSV.^2);
rsquaredOLS = 1 - RSSOLS/TSS;
yfittestOLSV = [ones(ntest,1) xtestM] * bOLSV; 
restestOLSV = ytestV-yfittestOLSV; 
RSStestOLS = sum(restestOLSV.^2);
rsquaredtestOLS = 1 - RSStestOLS/TSStest;
figure(1)
clf
plot(yV,yfitOLSV,'.')
hold on
plot(ytestV,yfittestOLSV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('OLS fit R^2=%1.4f and predict R^2=%1.4f',rsquaredOLS,rsquaredtestOLS))
figure(2)
clf
plot(yV,resOLSV/std(resOLSV),'.','Markersize',10)
hold on
plot(ytestV,restestOLSV/std(restestOLSV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('OLS, blue->fit, red->predict')

%% Stepwise fit
[bstepV,sdbV,pvalV,inmodel,stats]=stepwisefit(xM,yV);
bstep0 = stats.intercept;
bstepV = [bstep0; bstepV].*[1 inmodel]';
yfitstepV = [ones(n,1) xM] * bstepV; 
resstepV = yV-yfitstepV; 
RSSstep = sum(resstepV.^2);
rsquaredstep = 1 - RSSstep/TSS;
yfitteststepV = [ones(ntest,1) xtestM] * bstepV; 
resteststepV = ytestV-yfitteststepV; 
RSSteststep = sum(resteststepV.^2);
rsquaredteststep = 1 - RSSteststep/TSStest;
figure(3)
clf
plot(yV,yfitstepV,'.')
hold on
plot(ytestV,yfitteststepV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('stepwise fit R^2=%1.4f and predict R^2=%1.4f',rsquaredstep,rsquaredteststep))
figure(4)
clf
plot(yV,resstepV/std(resstepV),'.','Markersize',10)
hold on
plot(ytestV,resteststepV/std(resteststepV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('stepwise, blue->fit, red->predict')

%% PCR
lambdaV = zeros(r,1);
lambdaV(1:d) = 1;
bPCRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
bPCRV = [my - mxV*bPCRV; bPCRV];
yfitPCRV = [ones(n,1) xM] * bPCRV; 
resPCRV = yfitPCRV - yV;     % Calculate residuals
RSSPCR = sum(resPCRV.^2);
rsquaredPCR = 1 - RSSPCR/TSS;
yfittestPCRV = [ones(ntest,1) xtestM] * bPCRV; 
restestPCRV = ytestV-yfittestPCRV; 
RSStestPCR = sum(restestPCRV.^2);
rsquaredtestPCR = 1 - RSStestPCR/TSStest;
figure(5)
clf
plot(yV,yfitPCRV,'.')
hold on
plot(ytestV,yfittestPCRV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('PCR fit R^2=%1.4f and predict R^2=%1.4f',rsquaredPCR,rsquaredtestPCR))
figure(6)
clf
plot(yV,resPCRV/std(resPCRV),'.','Markersize',10)
hold on
plot(ytestV,restestPCRV/std(restestPCRV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('PCR, blue->fit, red->predict')

%% PLS
[Xloadings,Yloadings,Xscores,Yscores,bPLSV] = plsregress(xM,yV,d);
yfitPLSV = [ones(n,1) xM]*bPLSV;
resPLSV = yfitPLSV - yV;     % Calculate residuals
RSSPLS = sum(resPLSV.^2);
rsquaredPLS = 1 - RSSPLS/TSS;
yfittestPLSV = [ones(ntest,1) xtestM] * bPLSV; 
restestPLSV = ytestV-yfittestPLSV; 
RSStestPLS = sum(restestPLSV.^2);
rsquaredtestPLS = 1 - RSStestPLS/TSStest;
figure(7)
clf
plot(yV,yfitPLSV,'.')
hold on
plot(ytestV,yfittestPLSV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('PLS fit R^2=%1.4f and predict R^2=%1.4f',rsquaredPLS,rsquaredtestPLS))
figure(8)
clf
plot(yV,resPLSV/std(resPLSV),'.','Markersize',10)
hold on
plot(ytestV,restestPLSV/std(restestPLSV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('PLS, blue->fit, red->predict')

%% Ridge regression
mu = RSSOLS/(n-p);
sigmaV = diag(sigmaM);
lambdaV = sigmaV.^2 ./ (sigmaV.^2 + mu);
bRRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
bRRV = [my - mxV*bRRV; bRRV];
yfitRRV = [ones(n,1) xM] * bRRV; 
resRRV = yfitRRV - yV;     % Calculate residuals
RSSRR = sum(resRRV.^2);
rsquaredRR = 1 - RSSRR/TSS;
yfittestRRV = [ones(ntest,1) xtestM] * bRRV; 
restestRRV = ytestV-yfittestRRV; 
RSStestRR = sum(restestRRV.^2);
rsquaredtestRR = 1 - RSStestRR/TSStest;
figure(9)
clf
plot(yV,yfitRRV,'.')
hold on
plot(ytestV,yfittestRRV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('RR fit R^2=%1.4f and predict R^2=%1.4f',rsquaredRR,rsquaredtestRR))
figure(10)
clf
plot(yV,resRRV/std(resRRV),'.','Markersize',10)
hold on
plot(ytestV,restestRRV/std(restestRRV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('RR, blue->fit, red->predict')

%% LASSO 
[bM,fitinfo] = lasso(xcM,ycV);
lassoPlot(bM,fitinfo,'PlotType','Lambda','XScale','log');
lambda = input('Give lambda >');
[lmin, ilmin] = min(abs(fitinfo.Lambda - lambda));
bLASSOV = bM(:,ilmin);
bLASSOV = [my - mxV*bLASSOV; bLASSOV];
yfitLASSOV = [ones(n,1) xM] * bLASSOV; 
resLASSOV = yfitLASSOV - yV;     % Calculate residuals
RSSLASSO = sum(resLASSOV.^2);
rsquaredLASSO = 1 - RSSLASSO/TSS;
yfittestLASSOV = [ones(ntest,1) xtestM] * bLASSOV; 
restestLASSOV = ytestV-yfittestLASSOV; 
RSStestLASSO = sum(restestLASSOV.^2);
rsquaredtestLASSO = 1 - RSStestLASSO/TSStest;
figure(11)
clf
plot(yV,yfitLASSOV,'.')
hold on
plot(ytestV,yfittestLASSOV,'.r')
xlabel('y')
ylabel('$\hat{y}$','Interpreter','Latex')
title(sprintf('LASSO fit R^2=%1.4f and predict R^2=%1.4f',rsquaredLASSO,rsquaredtestLASSO))
figure(12)
clf
plot(yV,resLASSOV/std(resLASSOV),'.','Markersize',10)
hold on
plot(ytestV,restestLASSOV/std(restestLASSOV),'.r','Markersize',10)
plot(xlim,1.96*[1 1],'--c')
plot(xlim,-1.96*[1 1],'--c')
xlabel('y')
ylabel('e^*')
title('LASSO, blue->fit, red->predict')

fprintf('\t OLS \t step \t PCR \t PLS \t RR \t LASSO \n');
fprintf('\t Regression coefficient vectors \n');
disp([ bOLSV bstepV bPCRV bPLSV bRRV bLASSOV])
fprintf('\t fit R^2 \n');
disp([NaN rsquaredOLS rsquaredstep rsquaredPCR rsquaredPLS rsquaredRR rsquaredLASSO])
fprintf('\t predict R^2 \n');
disp([NaN rsquaredtestOLS rsquaredteststep rsquaredtestPCR rsquaredtestPLS rsquaredtestRR rsquaredtestLASSO])
