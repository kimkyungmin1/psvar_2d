% Apply instruments 
% IR for MP: Element at (1,2) is zero

function [mrot,mze] = fidinst(mb0,viv,mrv,mpv,nlag)

nrow = size(mrv)*[1;0];
nvar = size(mrv)*[0;1];

s_z = viv;
s_y = mrv((nlag + 1):nrow,:)-mpv;
mzy = s_z'*s_y;
mze = mzy/(mb0');

mrot = eye(nvar);
mrot(1:2, 3:nvar) = -mze(:,1:2) \ mze(:,3:nvar);
for ii = 3:nvar
    if ii > 3
        mrot(:, ii) = mrot(:, ii) - mrot(:, 3:(ii - 1)) * mrot(:, 3:(ii - 1))' * mrot(:, ii);
    end
    mrot(:, ii) = mrot(:, ii) / sqrt(mrot(:, ii)' * mrot(:, ii));
end

mrot(:, 1) = mrot(:, 1) - mrot(:, 3:nvar) * mrot(:, 3:nvar)' * mrot(:, 1);
mrot(:, 1) = mrot(:, 1) / sqrt(mrot(:, 1)' * mrot(:, 1));
mrot(:, 2) = mrot(:, 2) - mrot(:, [1 3:nvar]) * mrot(:, [1 3:nvar])' * mrot(:, 2);
mrot(:, 2) = mrot(:, 2) / sqrt(mrot(:, 2)' * mrot(:, 2));

mr2 = zeros(nvar, nvar);
mr2(1, 1) = 1;
mr2(2, 2) = 1;
mr2(3, 3) = 1;
cbb = mb0 * mrot;
for ii = 4:nvar
    mr2(ii, ii) = 1;
    cb = cbb * mr2;
    for ij = 3:(ii - 1)
        mr2(:, ii) = mr2(:, ii) - cb(ij, ii) / cb(ij, ij) * mr2(:, ij);
        cb = cbb * mr2;
    end
    mr2(:, ii) = mr2(:, ii) / sqrt(mr2(:, ii)' * mr2(:, ii));
end
for ii = 1:(nvar - 3)
    mr2(:, nvar - ii) = mr2(:, nvar - ii) - mr2(:, (nvar - ii + 1):nvar) * mr2(:, (nvar - ii + 1):nvar)' * mr2(:, nvar - ii);
    mr2(:, nvar - ii) = mr2(:, nvar - ii) / sqrt(mr2(:, nvar - ii)' * mr2(:, nvar - ii));
end

cb = cbb * mr2;

mr2(1, 2) = -cbb(1, 2) / cbb(1, 1);
mr2(:, 2) = mr2(:, 2) / sqrt(mr2(:, 2)' * mr2(:, 2));
mr2(:, 1) = mr2(:, 1) - mr2(:, 2) * mr2(:, 2)' * mr2(:, 1);
mr2(:, 1) = mr2(:, 1) / sqrt(mr2(:, 1)' * mr2(:, 1));

mrot = mrot * mr2;
cb = mb0 * mrot;
for ii = 1:nvar
    if cb(ii, ii) < 0
        mrot(:, ii) = -mrot(:, ii);
    end
end

end
    