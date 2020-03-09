%% example for Multivariate linear regression
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

tmpdata = table2array(data);
tmpdata(:,1) = [] % remove col with time
target = tmpdata(:,1); % set appliance as target
tmpdata(:,1) = [] % remove target col(appliances)
tmpdata(:,1) = [] % remove lights col
tmpdata = normalize(tmpdata) % necessary to normalize data before performing PCA. 
                             % all variables have the same standard deviation
%% Implement dimensionality reduction via pca 
numberOfDimensions = 3;
[coeff, score, latent, tsquared, explained] = pca(tmpdata);
reducedDimension = coeff(:,1:numberOfDimensions);
reducedData = tmpdata * reducedDimension;
tmpdata = reducedData;
%% Extract the response and predictor data
disp("ha");

regions = data.Properties.VariableNames;
regions = regions(:,1:numberOfDimensions);
Y = tmpdata
x = target
[n,d] = size(tmpdata);
disp('Hi');

%% Plot the data
figure;
plot(x,Y,'x')
legend(regions,'Location','NorthWest')

disp("haaaaaaaa");
%% Fit the multivariate regression model
X = cell(n,1);
for i = 1:n
	X{i} = [eye(d) repmat(x(i),d,1)];
end

% Sigma contains estimates of the -by- variance-covariance matrix , for the between-region concurrent correlations
% beta contains estimates of the -dimensional coefficient vector
[beta,Sigma] = mvregress(X,Y); % will throw error : 
                                % handle not positive covariance mvregress
% The covariance matrix is not positive definite because it is singular.
% That means that at least one of your variables can be expressed as a linear 
% combination of the others. You do not need all the variables as the value of 
% at least one can be determined from a subset of the others. I would suggest
% adding variables sequentially and checking the covariance matrix at each step. 
% If a new variable creates a singularity drop it and go on the the next one.
% Eventually you should have a subset of variables with a postive definite 
% covariance matrix.
% DONE
%% Plot the fitted regression model.

B = [beta(1:d)';repmat(beta(end),1,d)];
xx = linspace(.5,3.5)';
fits = [ones(size(xx)),xx]*B;

figure;
h = plot(x,Y,'x',xx,fits,'-');
for i = 1:d
	set(h(d+i),'color',get(h(i),'color'));
end
legend(regions,'Location','NorthWest');