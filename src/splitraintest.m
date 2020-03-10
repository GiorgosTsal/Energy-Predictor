%% clear env,get and set current directory
clc
clear
currdir = pwd
fprintf(currdir)
userpath(currdir) %set working directory to current dir of .m file
%% Sample data (54000 x 10)
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
tmpdata = table2array(data);
tmpdata(:,1) = []
tmpdata(:,end) = []
tmpdata(:,end) = []
[n, p] = size(tmpdata)

% Cross validation (train: 70%, test: 30%)
cv = cvpartition(size(tmpdata,1),'HoldOut',0.3);
idx = cv.test;

% Separate to training and test data
dataTrain = tmpdata(~idx,:);
dataTest  = tmpdata(idx,:);