%% clear env,get and set current directory
% run: close all to close all figures
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% load my dataset and convert date to number 
name = '/energydata_complete.csv';
filename = strcat(currdir,name)
data = importfile(filename)
data=data(1:6*24*63, :); % first 2,5 months (4.5months/2)
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
% tmpdata(:,1) = []; % remove first row
% tmpdata(:,end) = [];
% tmpdata(:,end) = [];
tmpdata = tmpdata(:,2:end); % remove first column(date)
tmpdata = tmpdata(:,1:(end-2)); % remove 2 last cols (random var 1 and random var 2)
%% split
[m,n] = size(tmpdata);
P = 0.70;
dataTrain = tmpdata(1:round(P*m),:);
dataTest = tmpdata((round(P*m)+1:end),:);
% idx = randperm(m)  ;
% dataTrain = tmpdata(idx(1:round(P*m)),:) ; 
% dataTest = tmpdata(idx(round(P*m)+1:end),:) ;

disp('john');
%% set target
yV=dataTrain(:,1); % target energy from appliances
xM = dataTrain(:,2:end); % without first column target(appliances)

%% run regressions
  alpha = 0.05;
    zalpha = norminv(1-alpha/2);
    [n,m] = size(xM);

    x=xM;
    y=yV;
    
%     x = zeros(n,m);
%     for i=1:m
%         x(:,i) =fitAR(xM(:,i),2);
%     end
%     y = fitAR(yV,2);
%     
    ytrain = y(1:floor(0.7*n));
    ytest = y(ceil(0.7*n):n);
    xtrain = x(1:floor(0.7*n),:);
    xtest = x(ceil(0.7*n):n,:);
    p=m;
    [n,m]=size(xtrain);
    TSS = sum((ytrain-mean(ytrain)).^2);
    mxV = mean(xtrain);
    xcM = xtrain - repmat(mxV,n,1); % centered data matrix
    my = mean(ytrain);
    ycV = ytrain - my;
    
    [uM,sigmaM,vM] = svd(xcM,'econ');
    r = size(sigmaM,1);
   
    %% OLS
    bOLSV = vM * inv(sigmaM) * uM'* ycV;
    % yfitOLSV = xcM * bOLSV + my; 
    bOLSV = [my - mxV*bOLSV; bOLSV];
    yfitOLSV = [ones(n,1) xtrain] * bOLSV; 
    resOLSV = ytrain-yfitOLSV; 
    RSSOLS = sum(resOLSV.^2);
    mu = RSSOLS/(n-p);
    rsquaredOLS = 1 - RSSOLS/TSS;
    figure(1)
    clf
    plot(ytrain,yfitOLSV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('OLS R^2=%1.4f',rsquaredOLS))
    figure(2)
    clf
    plot(ytrain,resOLSV/std(resOLSV),'.','Markersize',10)
    hold on
    plot(xlim,zalpha*[1 1],'--c')
    plot(xlim,-zalpha*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('OLS')
    %test    
    [n,m]=size(xtest);
    TSS = sum((ytest-mean(ytest)).^2);

    yfitOLSV = [ones(n,1) xtest] * bOLSV; 
    resOLSV = ytest-yfitOLSV; 
    RSSOLS = sum(resOLSV.^2);
    rsquaredOLS = 1 - RSSOLS/TSS;
    figure(3)
    clf
    plot(ytest,yfitOLSV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('OLS R^2=%1.4f',rsquaredOLS))
    figure(4)
    clf
    plot(ytest,resOLSV/std(resOLSV),'.','Markersize',10)
    hold on
    plot(xlim,zalpha*[1 1],'--c')
    plot(xlim,-zalpha*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('OLS')
    
    %% PCR
    
    % dimension reduction
    d=m;
    [n,m]=size(xtrain);
    lambdaV = zeros(r,1);
    lambdaV(1:d) = 1;
    bPCRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
    bPCRV = [my - mxV*bPCRV; bPCRV];
    yfitPCRV = [ones(n,1) xtrain] * bPCRV; 
    resPCRV = yfitPCRV - ytrain;     % Calculate residuals
    RSSPCR = sum(resPCRV.^2);
    rsquaredPCR = 1 - RSSPCR/TSS;
    figure(5)
    clf
    plot(ytrain,yfitPCRV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('PCR R^2=%1.4f',rsquaredPCR))
    figure(6)
    clf
    plot(ytrain,resPCRV/std(resPCRV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('PCR')
    
    %test
    [n,m]=size(xtest);
    yfitPCRV = [ones(n,1) xtest] * bPCRV; 
    resPCRV = yfitPCRV - ytest;     % Calculate residuals
    RSSPCR = sum(resPCRV.^2);
    rsquaredPCR = 1 - RSSPCR/TSS;
    figure(7)
    clf
    plot(ytest,yfitPCRV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('PCR R^2=%1.4f',rsquaredPCR))
    figure(8)
    clf
    plot(ytest,resPCRV/std(resPCRV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('PCR')
    
    
    %% Stepwise fit
    [n,m]=size(xtrain);
    
    [bstepV,sdbV,pvalV,inmodel,stats]=stepwisefit(xtrain,ytrain);
    bstep0 = stats.intercept;
    bstepV = [bstep0; bstepV].*[1 inmodel]';
    yfitstepV = [ones(n,1) xtrain] * bstepV; 
    resstepV = ytrain-yfitstepV; 
    RSSstep = sum(resstepV.^2);
    rsquaredstep = 1 - RSSstep/TSS;
    figure(9)
    clf
    plot(ytrain,yfitstepV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('stepwise regression R^2=%1.4f',rsquaredstep))
    figure(10)
    clf
    plot(ytrain,resstepV/std(resstepV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('stepwise regression')

    %test 
    
    [n,m]=size(xtest);
    
    yfitstepV = [ones(n,1) xtest] * bstepV; 
    resstepV = ytest-yfitstepV; 
    RSSstep = sum(resstepV.^2);
    rsquaredstep = 1 - RSSstep/TSS;
    figure(11)
    clf
    plot(ytest,yfitstepV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('stepwise regression R^2=%1.4f',rsquaredstep))
    figure(12)
    clf
    plot(ytest,resstepV/std(resstepV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('stepwise regression')
    
    
    %% PLS
    [n,m]=size(xtrain);
    [Xloadings,Yloadings,Xscores,Yscores,bPLSV] = plsregress(xtrain,ytrain,d);
    yfitPLSV = [ones(n,1) xtrain]*bPLSV;
    resPLSV = yfitPLSV - ytrain;     % Calculate residuals
    RSSPLS = sum(resPLSV.^2);
    rsquaredPLS = 1 - RSSPLS/TSS;
    figure(13)
    clf
    plot(ytrain,yfitPLSV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('PLS R^2=%1.4f',rsquaredPLS))
    figure(14)
    clf
    plot(ytrain,resPLSV/std(resPLSV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('PLS')
    
    %test 
    [n,m]=size(xtest);
   
    yfitPLSV = [ones(n,1) xtest]*bPLSV;
    resPLSV = yfitPLSV - ytest;     % Calculate residuals
    RSSPLS = sum(resPLSV.^2);
    rsquaredPLS = 1 - RSSPLS/TSS;
    figure(15)
    clf
    plot(ytest,yfitPLSV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('PLS R^2=%1.4f',rsquaredPLS))
    figure(16)
    clf
    plot(ytest,resPLSV/std(resPLSV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('PLS')
    
    
    %% Ridge
    [n,m]=size(xtrain);
    p=m;
    sigmaV = diag(sigmaM);
    lambdaV = sigmaV.^2 ./ (sigmaV.^2 + mu);
    bRRV = vM * diag(lambdaV) * inv(sigmaM) * uM'* ycV;
    bRRV = [my - mxV*bRRV; bRRV];
    %bRRV = ridge(ytrain,xtrain,1e-3);
    yfitRRV = [ones(n,1) xtrain]* bRRV; 
    resRRV = yfitRRV - ytrain;     % Calculate residuals
    RSSRR = sum(resRRV.^2);
    rsquaredRR = 1 - RSSRR/TSS;
    figure(17)
    clf
    plot(ytrain,yfitRRV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('RR R^2=%1.4f',rsquaredRR))
    figure(18)
    clf
    plot(ytrain,resRRV/std(resRRV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('RR')

    
    % test 
    
    [n,m]=size(xtest);
    
    yfitRRV = [ones(n,1) xtest] * bRRV; 
    resRRV = yfitRRV - ytest;     % Calculate residuals
    RSSRR = sum(resRRV.^2);
    rsquaredRR = 1 - RSSRR/TSS;
    figure(19)
    clf
    plot(ytest,yfitRRV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('RR R^2=%1.4f',rsquaredRR))
    figure(20)
    clf
    plot(ytest,resRRV/std(resRRV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('RR')

    %% Lasso
    [n,m]=size(xtrain);
    [bM,fitinfo] = lasso(xcM,ycV);
    lassoPlot(bM,fitinfo,'PlotType','Lambda','XScale','log');
    lambda = input('Give lambda >');
    [lmin, ilmin] = min(abs(fitinfo.Lambda - lambda));
    bLASSOV = bM(:,ilmin);
    bLASSOV = [my - mxV*bLASSOV; bLASSOV];
    yfitLASSOV = [ones(n,1) xtrain] * bLASSOV; 
    resLASSOV = yfitLASSOV - ytrain;     % Calculate residuals
    RSSLASSO = sum(resLASSOV.^2);
    rsquaredLASSO = 1 - RSSLASSO/TSS;
    figure(22)
    clf
    plot(ytrain,yfitLASSOV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('LASSO R^2=%1.4f',rsquaredLASSO))
    figure(23)
    clf
    plot(ytrain,resLASSOV/std(resLASSOV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('LASSO')
    
    % test
    
    [n,m]=size(xtest);
    
    yfitLASSOV = [ones(n,1) xtest] * bLASSOV; 
    resLASSOV = yfitLASSOV - ytest;     % Calculate residuals
    RSSLASSO = sum(resLASSOV.^2);
    rsquaredLASSO = 1 - RSSLASSO/TSS;
    figure(24)
    clf
    plot(ytest,yfitLASSOV,'.')
    hold on
    xlabel('y')
    ylabel('$\hat{y}$','Interpreter','Latex')
    title(sprintf('LASSO R^2=%1.4f',rsquaredLASSO))
    figure(25)
    clf
    plot(ytest,resLASSOV/std(resLASSOV),'.','Markersize',10)
    hold on
    plot(xlim,1.96*[1 1],'--c')
    plot(xlim,-1.96*[1 1],'--c')
    xlabel('y')
    ylabel('e^*')
    title('LASSO')
    
    fprintf('\t OLS \t step \t PCR \t PLS \t RR \t LASSO \n');
    fprintf('\t Regression coefficient vectors \n');
    disp([bOLSV bstepV bPCRV bPLSV bRRV bLASSOV])
    fprintf('\t R^2 values \n');
    disp([rsquaredstep rsquaredOLS rsquaredPCR rsquaredPLS rsquaredRR rsquaredLASSO])