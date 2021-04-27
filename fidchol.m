% Compute structural shocks 
function mb0 = fidchol(mrv,mpv,nlag)

nrow = size(mrv)*[1;0];
me0 = mrv((nlag+1):nrow,:)-mpv;
sigma0 = me0'*me0 / (nrow-nlag);
mb0 = chol(sigma0, 'lower');

end
    