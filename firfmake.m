% Compute impulse responses 
function irf0 = firfmake(mb0,mcf1,nmaxt)

mcf1_s = size(mcf1);
mcf1_sd = size(mcf1_s)*[0;1];
if mcf1_sd == 2
    nlag=1;
else
    nlag=mcf1_s(3);
end
nvar = size(mb0)*[1;0];
irf0 = zeros(nmaxt+1,nvar,nvar);
rx0 = zeros(nlag+nmaxt+1,nvar);
for ik = 1:nvar
    rx0(nlag+1,:) = mb0(:,ik);
    for ii = 1:nmaxt
        ry = zeros(nvar,1);
        for ij = 1:nlag
            rb = mcf1(:,:,ij);
            ry = ry + rb * rx0(nlag+ii-ij+1,:)';
        end
        rx0(nlag+ii+1,:) = ry;
    end
    irf0(:, :, ik) = rx0((nlag + 1):(nlag + nmaxt + 1), :);
end

end
    