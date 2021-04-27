% Compute VAR coefficient matrices

function [mcf1,mcon1,rx] = fvarcoef(vars,nlag)

nrow = size(vars)*[1;0];
nvar = size(vars)*[0;1];
rxvars = zeros(nrow-nlag,nlag*nvar);
for ii = 1:nvar
    for ij = 1:nlag
        rxvars(:,(ii-1)*nlag+ij) = vars((nlag+1-ij):(nrow-ij),ii);
    end
end
mcf1 = zeros(nvar,nvar,nlag);
mcon1 = zeros(nvar,1);
rx = [ones(nrow-nlag,1) rxvars];

for ii = 1:nvar
    ry = vars((nlag + 1):nrow,ii);
    rb = rx\ry;
    mcon1(ii) = rb(1);
    for ij = 1:nlag
        mcf1(ii,:,ij) = rb((1+ij):nlag:(1+nlag*(nvar-1)+ij));
    end
end

end
    