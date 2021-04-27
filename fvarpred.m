% Compute predicted values for the VAR

function mpv = fvarpred(mcf,mcon,rxall)

nrow = size(rxall)*[1;0];
smcf = size(mcf);
if (size(smcf)*[0;1]) == 2
    nvar = size(mcf)*[1;0];
    nlag = 1;
else
    nvar = size(mcf)*[1;0;0];
    nlag = size(mcf)*[0;0;1];
end
mpv = zeros(nrow,nvar);
for ii = 1:nvar
    rb = zeros(1+nvar*nlag,1);
    rb(1) = mcon(ii);
    for ij = 1:nlag
        rb((1+ij):nlag:(1+(nvar-1)*nlag+ij)) = mcf(ii,:,ij);
    end
    mpv(:,ii) = rxall*rb;
end

end
    