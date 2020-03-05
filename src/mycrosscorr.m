function ccV = mycrosscorr(xV,yV,tau)
% ccV = mycrosscorr(xV,yV,tau)
% Calls the function crosscorr when 'tau>0' and the function corrcoef when 
% 'tau=0'.

if tau==0
    tmpM = corrcoef(xV,yV);
    ccV = tmpM(1,2);
else
    ccV = crosscorr(xV,yV,tau);
end
