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

data2.date = datenum(data2.date, 'yyyy-mm-dd HH:MM:SS');
ts = data2.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data2.date = ts;

tmpdata = table2array(data2);
tmpdata = tmpdata(:,2:end); % remove first column(date)
%% Plot the distribution of features
%assign variable names
nameM = data2.Properties.VariableNames;
nameM = nameM';
[n,m] = size(tmpdata);
for i=1:m
    makis = cell2mat(nameM(i,:));
    title(makis)  
    h = figure(i);
    histogram(tmpdata(:,i))
    saveas(h,makis,'png')
end
%% plot correlation between all variables with each other(lags AF)
%corrplot(data)