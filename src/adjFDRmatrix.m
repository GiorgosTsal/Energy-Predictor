function adjM = adjFDRmatrix(xM,alpha,symmetric)
% adjM = adjFDRmatrix(xM,alpha,symmetric)
% Function adjFDRmatrix creates an adjacency matrix of zeros and ones from
% a given square KxK matrix 'xM' of p-values, where ones regards
% significant p-values according to the False Discovery Rate (FDR)
% criterion for a given 'alpha'. If symmetric==1, then 'xM' is supposed to
% be symmetric and only the upper (or lower) triangular off diagonal matrix
% is considered (K*(K-1)/2 values), otherwise all but the diagonal
% components are considered (K*(K-1) values).
% INPUT 
% - xM          : square K x K matrix of p-values.
% - alpha       : the significance limit for FDR
% - symmetric   : if 1, xM is to be considered as symmetric (default ~=1)
% OUTPUT
% - adjM        : the adjacency K x K matrix zeros and ones (ones for
%                 significant rejection, small p-values). 

if nargin==2
    symmetric = 0;
elseif nargin ==1
    symmetric = 0;
    alpha = 0.05;
end
if isempty(symmetric), symmetric=0; end
if isempty(alpha), alpha=0.05; end

[K,K1] = size(xM);
if K~=K1
    error('The input matrix of p-values must be square.');
end

if symmetric
    m = K*(K-1)/2;
    xvecM = NaN*ones(m,3);
    count = 0;
    for i=1:K
        jV = [i+1:K]';
        nj = length(jV);
        xvecM(count+[1:nj]',1)= i*ones(nj,1);
        xvecM(count+[1:nj]',2)= jV;
        xvecM(count+[1:nj]',3)= xM(i,jV)';
        count = count+nj;
    end
else
    m = K*(K-1);
    xvecM = NaN*ones(m,3);
    count = 0;
    for i=1:K
        jV = [1:i-1 i+1:K]';
        nj = length(jV);
        xvecM(count+[1:nj]',1)= i*ones(nj,1);
        xvecM(count+[1:nj]',2)= jV;
        xvecM(count+[1:nj]',3)= xM(i,jV)';
        count = count+nj;
    end
end
[oxvecV, ixvecV]=sort(xvecM(:,3));
iV = find(oxvecV<=alpha*[1:m]'/m);
adjM = zeros(K,K);
if ~isempty(iV)
    ixvecV = ixvecV(1:iV(end));
    for i=1:length(ixvecV)
        adjM(xvecM(ixvecV(i),1),xvecM(ixvecV(i),2))=1;
        if symmetric
            adjM(xvecM(ixvecV(i),2),xvecM(ixvecV(i),1))=1;
        end    
    end
end
