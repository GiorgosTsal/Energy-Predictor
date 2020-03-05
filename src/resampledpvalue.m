function pval = resampledpvalue(xV,twosided)
% pval = resampledpvalue(xV,twosided)
% Computes the p-value for a resampling test, where M+1 statistic values
% are given in 'xV', the first is the statistic value for the original data
% set, and the rest M+1 for the surrogates. If the test is two-sided then
% twosided=1 (default), otherwise left or right sided test is applied.
% INPUTS 
% - xV          : the array of M+1 statistic values (first for the original data set)
% - twosided    : if 1-> twosided test (default), if -1-> left sided,
%                 otherwise right sided test. 
% OUTPUT 
% - pval        : the p-value of the resampling test.

if nargin==1
    twosided = 1;
end

nx = length(xV);
[oxV,ixV]=sort(xV);
rnkx = find(ixV == 1);
isamexV = find(xV==xV(1));
if length(isamexV)==nx
    rnkx=round(nx/2);
elseif length(isamexV)>=2
    irand = ceil(length(isamexV)*rand);
    rnkx = rnkx+irand-1;
end
if twosided==1
    if rnkx > 0.5*nx
        pval = 2*(1-(rnkx-0.326)/(nx+0.348));
    else
        pval = 2*(rnkx-0.326)/(nx+0.348);
    end
elseif twosided == -1
    pval = (rnkx-0.326)/(nx+0.348);
else
    pval = (1-(rnkx-0.326)/(nx+0.348));
end

