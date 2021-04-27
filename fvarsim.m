% Construct counterfactual paths
function mcy = fvarsim(vars,mcon1,mcf1,nlag,pmu)

nrow = size(vars)*[1;0];
nvar = size(vars)*[0;1];
mcy = zeros(nrow,nvar);
mcy(1:nlag,:) = vars(1:nlag,:);

for ii = (nlag+1):nrow
    ry = mcon1+pmu(ii-nlag,:)';
    for ij = 1:nlag
        ry = ry + mcf1(:,:,ij)*mcy(ii-ij,:)';
    end
    mcy(ii,:) = ry;
end

end
    