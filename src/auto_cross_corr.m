%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% Lab 3: Explore autocorrelations and cross correlations in greek stocks
maxtau = 20;
maxtau2 = 10;
p = 1;
M = 100;
alpha = 0.05;
zalpha = 1.96;

rng(1);
%% load my dataset and convert date to number 
name = '/energydata_complete.csv';
filename = strcat(currdir,name)
data = importfile(filename)

data(:,1) = []; %remove col with timestamps
data(:,end) = []; %remove ranodm var1
data(:,end) = []; %remove random var2

% data.date = datenum(data.date, 'yyyy-mm-dd HH:MM:SS');
% ts = data.date; % temp variable 
% ts = ts*24*60*60; % tranform date to seconds
% ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
% data.date = ts;
disp('Hi');

%assign variable names
nameM = data.Properties.VariableNames;
nameM = nameM';
tmpdata = table2array(data)

[n,m]=size(tmpdata);
yM = tmpdata';
%% First index
indx1 = ceil(m*rand);
%y1V = yM(:,indx1);
y1V = yM(1,:)
y1V = y1V'
name1 = cell2mat(nameM(indx1,:));
% If NaN replace them with interpolated values
i1V = find(isnan(y1V));
if ~isempty(i1V)
    iokV = setdiff([1:n]',i1V);
    y1V(i1V) = interp1(iokV,y1V(iokV),i1V,'spline');
end

figure(1)
clf
plot(y1V,'.-')
xlabel('day t')
ylabel('y(t)')
title(sprintf('variable %s',name1))

acy1M = autocorrelation(y1V,maxtau);
figure(2)
clf
plot(acy1M(:,1),acy1M(:,2),'.-')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s',name1))

%% Second index
indx2 = ceil(m*rand);
%y2V = yM(:,indx2);
y2V = yM(1,:)
y2V = y2V'
name2 = cell2mat(nameM(indx2,:));
% If NaN replace them with interpolated values
i2V = find(isnan(y2V));
if ~isempty(i2V)
    iokV = setdiff([1:n]',i2V);
    y2V(i2V) = interp1(iokV,y2V(iokV),i2V,'spline');
end

figure(3)
clf
plot(y2V,'.-')
xlabel('day t')
ylabel('y(t)')
title(sprintf('variable %s',name2))

acy2M = autocorrelation(y2V,maxtau);
figure(4)
clf
plot(acy2M(:,1),acy2M(:,2),'.-')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s',name2))

ccyV = mycrosscorr(y1V,y2V,maxtau2);
figure(5)
clf
plot([-maxtau2:maxtau2]',ccyV,'.-')
xlabel('lag \tau')
ylabel('r_{XY}(\tau)')
title(sprintf('Cross Corr of variable %s and %s',name1,name2))

%% Prewhitening, first method: Fit AR model and use the residuals
ey1V = fitAR(y1V,p);
figure(6)
clf
plot(ey1V,'.-')
xlabel('day t')
ylabel('residual of AR of y(t)')
title(sprintf('variable %s, residual from AR(%d)',name1,p))

acey1M = autocorrelation(ey1V,maxtau);
figure(7)
clf
plot(acey1M(:,1),acey1M(:,2),'.-')
hold on
plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s, residual from AR(%d)',name1,p))

ey2V = fitAR(y2V,p);
figure(8)
clf
plot(ey2V,'.-')
xlabel('day t')
ylabel('residual of AR of y(t)')
title(sprintf('variable %s, residual from AR(%d)',name2,p))

acey2M = autocorrelation(ey2V,maxtau);
figure(9)
clf
plot(acey2M(:,1),acey2M(:,2),'.-')
hold on
plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s, residual from AR(%d)',name2,p))

cceyV = mycrosscorr(ey1V,ey2V,maxtau2);
figure(10)
clf
plot([-maxtau2:maxtau2]',cceyV,'.-')
hold on
plot([-maxtau2 maxtau2],(zalpha/sqrt(n))*[1 1],'c--')
plot([-maxtau2 maxtau2],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_{XY}(\tau)')
title(sprintf('Cross Corr of variable %s and %s, residual from AR(%d)',name1,name2,p))

%% Prewhitening, second method: Take first differences (or log returns)
ey1V = log(y1V(2:n))-log(y1V(1:n-1));
figure(11)
clf
plot(ey1V,'.-')
xlabel('day t')
ylabel('log returns of y(t)')
title(sprintf('variable %s, log returns',name1))

acey1M = autocorrelation(ey1V,maxtau);
figure(12)
clf
plot(acey1M(:,1),acey1M(:,2),'.-')
hold on
plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s, log returns',name1))

ey2V = log(y2V(2:n))-log(y2V(1:n-1));
figure(13)
clf
plot(ey2V,'.-')
xlabel('day t')
ylabel('log returns of y(t)')
title(sprintf('variable %s, log returns',name2))

acey2M = autocorrelation(ey2V,maxtau);
figure(14)
clf
plot(acey2M(:,1),acey2M(:,2),'.-')
hold on
plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of variable %s, log returns',name2))

cceyV = mycrosscorr(ey1V,ey2V,maxtau2);
figure(15)
clf
plot([-maxtau2:maxtau2]',cceyV,'.-')
hold on
plot([-maxtau2 maxtau2],(zalpha/sqrt(n))*[1 1],'c--')
plot([-maxtau2 maxtau2],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_{XY}(\tau)')
title(sprintf('Cross Corr of variable %s and %s, log returns',name1,name2))

