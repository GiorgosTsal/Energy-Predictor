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
data1=data(1:6*24*63, :)
data2=data(6*24*63+1:6*24*63*2,:)

data1.date = datenum(data1.date, 'yyyy-mm-dd HH:MM:SS');
ts = data1.date; % temp variable 
ts = ts*24*60*60; % tranform date to seconds
ts = ts - ts(1); % subtract sample one from all the other time samples(to start from zero secs)
data1.date = ts;

%assign variable names
nameM = data1.Properties.VariableNames;
%%nameM = nameM';
tmpdata = table2array(data1)

[n,m]=size(tmpdata);
%yM = tmpdata';
yM=tmpdata;

%{
%% First index
for indx1 = 1:m

    %y1V = yM(:,indx1);
    y1V = yM(:,indx1);

    
    name1 = cell2mat(nameM(indx1));
    % If NaN replace them with interpolated values
    i1V = find(isnan(y1V));
    if ~isempty(i1V)
        iokV = setdiff([1:n]',i1V);
        y1V(i1V) = interp1(iokV,y1V(iokV),i1V,'spline');
    end

    figure(2*indx1-1)
    clf
    plot(y1V,'.-')
    xlabel('day t')
    ylabel('y(t)')
    title(sprintf('variable %s',name1))

    acy1M = autocorrelation(y1V,maxtau);
    figure(2*indx1)
    clf
    plot(acy1M(:,1),acy1M(:,2),'.-')
    xlabel('lag \tau')
    ylabel('r_Y(\tau)')
    title(sprintf('autocorrelation of variable %s',name1))
end


%}
%%crosscorrelation with noise


%{
for indx1=1:m-1
    for indx2=indx1+1:m
        y1V = yM(:,indx1);
        y2V = yM(:,indx2);
        name1 = cell2mat(nameM(indx1));
        name2 = cell2mat(nameM(indx2));
        ccyV = mycrosscorr(y1V,y2V,maxtau2);
        figure(5)
        clf
        plot([-maxtau2:maxtau2]',ccyV,'.-')
        xlabel('lag \tau')
        ylabel('r_{XY}(\tau)')
        title(sprintf('Cross Corr of stock %s and %s',name1,name2))

    end
end
%}



%%autocorrelation after whittening with method 1
%{
for indx1=1:m
    y1V = yM(:,indx1);
    ey1V = fitAR(y1V,p);
    figure(6)
    clf
    plot(ey1V,'.-')
    xlabel('day t')
    ylabel('residual of AR of y(t)')
    name1 = cell2mat(nameM(indx1));
    title(sprintf('Variable %s, residual from AR(%d)',name1,p))

    acey1M = autocorrelation(ey1V,maxtau);
    figure(7)
    clf
    plot(acey1M(:,1),acey1M(:,2),'.-')
    hold on
    plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
    plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
    xlabel('lag \tau')
    ylabel('r_Y(\tau)')
    title(sprintf('autocorrelation of stock %s, residual from AR(%d)',name1,p))
end
%}




%%cross correlation after whittening with method 1
%{
for indx1=1:m-1
    for indx2=indx1+1:m
        y1V = yM(:,indx1);
        y2V = yM(:,indx2);
        ey1V = fitAR(y1V,p);
        ey2V = fitAR(y2V,p);
        name1 = cell2mat(nameM(indx1));
        name2 = cell2mat(nameM(indx2));
        
        
        cceyV = mycrosscorr(ey1V,ey2V,maxtau2);
        figure(10)
        clf
        plot([-maxtau2:maxtau2]',cceyV,'.-')
        hold on
        plot([-maxtau2 maxtau2],(zalpha/sqrt(n))*[1 1],'c--')
        plot([-maxtau2 maxtau2],-(zalpha/sqrt(n))*[1 1],'c--')
        xlabel('lag \tau')
        ylabel('r_{XY}(\tau)')
        title(sprintf('Cross Corr of variables %s and %s, residual from AR(%d)',name1,name2,p))

    end
end
%}


%%autocorrelation after whittening with method 2
%{
for indx1=1:m
    y1V = yM(:,indx1);
    ey1V = log(y1V(2:n))-log(y1V(1:n-1));

    name1 = cell2mat(nameM(indx1));
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
end
%}

%%cross correlation after whittening with method 1
i=1
for indx1=1:m-1
    for indx2=indx1+1:m
        y1V = yM(:,indx1);
        y2V = yM(:,indx2);
        ey1V = log(y1V(2:n))-log(y1V(1:n-1));
        ey2V = log(y2V(2:n))-log(y2V(1:n-1));
        name1 = cell2mat(nameM(indx1));
        name2 = cell2mat(nameM(indx2));
        
        
        cceyV = mycrosscorr(ey1V,ey2V,maxtau2);
        figure(i)
        clf
        plot([-maxtau2:maxtau2]',cceyV,'.-')
        hold on
        plot([-maxtau2 maxtau2],(zalpha/sqrt(n))*[1 1],'c--')
        plot([-maxtau2 maxtau2],-(zalpha/sqrt(n))*[1 1],'c--')
        xlabel('lag \tau')
        ylabel('r_{XY}(\tau)')
        title(sprintf('Cross Corr of stock %s and %s, log returns',name1,name2))
        i=i+1
    end
end


