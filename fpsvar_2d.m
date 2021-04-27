% PURPOSE: Run SVAR with instruments.

% INSTRUMENTS: 2 shocks, 2 instruments.
% (IF YOU HAVE JUST 1 INSTRUMENT, PUT ANY SERIES AS THE SECOND INSTRUMENT
% AND SET IDCOND=2. YOU STILL HAVE THE CORRECT RESULT FOR THE FIRST INSTRUMENT)

% IDENTIFYING ASSUMPTIONS: 
% IF IDCOND not equal 1, ~(IDCOND==1):
% 1. By design, VAR with n variables identifies n structural shocks.
% 2. Only the first structural shock has non-zero covariance with the first
% instrument.
% 3. Only the first and second structural shocks have non-zero covariance
% with the first two instruments.
% 4. In most applications, you do not need to worry about the other shocks.
% They are identifiend as follows: shock n has zero t=0 impact on shocks 3, 4,
% ..., n-1; shock n-1 has zero t=0 impact on shocks 3, 4, ..., n-2; and so on.
% IF IDCOND == 1:
% Change in step 2 only: The t=0 impact of the second structural shock on the
% first variable is zero.

% IRF scaling: Shock n's t=0 impact on variable n is 1. To
% scale it so that shock n's t=0 impact on variable n is a, you can
% multiply IRF(:,n,:), CFBD{1}(:,n,:), and CFBD{2}(:,n,:) by a.

% INPUT:
% VARS: T*N matrix of VAR variables
% INST: T*2 matrix of instruments.
% (1. The program does not demean the instrument. You need to enter demeaned
% instruments. 2. You can set missing instruments as zero. 3. If you want to 
% demean the instruments in each simulation, you can easily do that. However, 
% you need to set the values of missing observations to be zero in each
% simulated path. It is not enough to set them to be zero in the input because
% then their values will no longer be zero due to demeaning.)
% NLAG: Number of lags in VAR.
% MAXHOR: Maximum horizon to calculate IRFs.
% NSIM: Number of bootstrap simulations to run.
% PCTL: A column vector of percentiles for confidence bands. For example,
% [0.05;0.95] will calculated 5th and 95th percentiles.
% IDCOND: If 1 identify with zero restriction. Otherwise, identify using
% restriction on the convariance between instruments and shocks.
% BLLEN: Bootstrap block length. Setting this to 1 results in sampling 
% with replacement. Setting this to a number greater than 1 results in 
% a moving block bootstrap as in Jentsch and Lunsford (2019) (except that 
% this code does not demean. In my experience this made hardly any
% difference, but this may depend on data).

% OUTPUT:
% IRF: (MAXHOR+1)*N*N matrix of IRFs. IRF(:,:,1) is the IRF to the first
% shock. The first dimension covers the time horizoon of IRF. (index 1: t=0,
% index 2: t=1, ..., index MAXHOR+1: t=MAXHOR). The second dimension covers
% the N variables in VAR. The third dimension covers the N structural 
% shocks.
% CFBD: Confidence bands. CFBD{ii} has the same shape as the IRF and
% represents pointwise PCTL(ii) percentiles.
% MISCOUT: other output variables of interest.
% MISCOUT.SSHOCK: T*N matrix of structural shocks. The first NLAG shocks
% are set to zero.
% MISCOUT.IZE: A 2*N matrix representing E[INST'*STRUCTURAL SHOCKS].


function [IRF,CFBD,MISCOUT] = fpsvar_2d(VARS,INST,NLAG,MAXHOR,NSIM,PCTL,IDCOND,BLLEN)

nrow = size(VARS)*[1;0];
nvar = size(VARS)*[0;1];

[mcf1,mcon1,rxall] = fvarcoef(VARS,NLAG); % Estimate VAR coefficients
mpv = fvarpred(mcf1,mcon1,rxall); % Calculated VAR predicted values (as OLS)
mu = VARS((NLAG+1):nrow,:)-mpv;

pmu = cell(NSIM+1,1);
piv = cell(NSIM+1,1);
zb = INST((NLAG+1):nrow,:);

for ii = 1:NSIM % Generate simulated paths for residuals and instruments
    bln = ceil((nrow-NLAG)/BLLEN);
    jk = randsample(nrow-NLAG-BLLEN+1,bln,true);
    ji = zeros(BLLEN*bln,1);
    for ij = 1:bln
        ji(((ij-1)*BLLEN+1):(ij*BLLEN)) = jk(ij):(jk(ij)+BLLEN-1);
    end
    ji = ji(1:(nrow-NLAG));
    pmu{ii} = mu(ji,:);
    piv{ii} = zb(ji,:);    
end

pmu{NSIM+1} = mu; % Save the original
piv{NSIM+1} = zb; % Save the original

pirf0 = cell(NSIM+1,1);
% Change for to parfor for parallel
for ii = 1:(NSIM+1) % Run VAR for each simulated path
    crv = fvarsim(VARS,mcon1,mcf1,NLAG,pmu{ii});
    [ccf,ccon,cxall] = fvarcoef(crv,NLAG);
    cpv = fvarpred(ccf,ccon,cxall);
    
    cb0 = fidchol(crv,cpv,NLAG); % Cholesky decomposition of the covariance matrix
    [crot,cze] = fidinst(cb0,piv{ii},crv,cpv,NLAG); % Apply instruments
        
    % Apply a further rotation for identification if ~(IDCOND==1)
    if IDCOND==1
        cb0iv = cb0*crot;
        ize = cze*crot;
    else
        cb0iv = cb0*crot;
        ize = cze*crot;
        crotdn = sqrt(ize(1,1)^2+ize(1,2)^2);
        crot2 = zeros(2,2);
        crot2(1,2) = ize(1,2)/crotdn;
        crot2(2,2) = -ize(1,1)/crotdn;
        crot2(1,1) = ize(1,1)/crotdn;
        crot2(2,1) = ize(1,2)/crotdn;
        crot2big = eye(nvar);
        crot2big(1:2,1:2) = crot2;
        cb0iv = cb0iv*crot2big;
        ize = ize*crot2big;
        if cb0iv(1,1) < 0
            cb0iv(:,1) = -cb0iv(:,1);
            ize(:,1) = -ize(:,1);
        end
        if cb0iv(2,2) < 0
            cb0iv(:,2) = -cb0iv(:,2);
            ize(:,2) = -ize(:,2);
        end
    end
    
    cirf0 = firfmake(cb0iv,ccf,MAXHOR);
    pirf0{ii} = cirf0;
    
    if ii==(NSIM + 1) % Create MISCOUT
        dcres = zeros(nrow,nvar);
        dcres((NLAG+1):nrow,:) = crv((NLAG+1):nrow,:)-cpv;
        MISCOUT.SSHOCK = dcres / (cb0iv');
        MISCOUT.ZE = ize;        
    end
end

npct = size(PCTL)*[1;0];
CFBD = cell(npct,1);
for ii = 1:npct
    crb = zeros(MAXHOR+1,nvar,nvar);
    for ij = 1:nvar
        for ik = 1:nvar
            for il = 1:(MAXHOR+1)
                pointall = zeros(NSIM,1);
                for ii2 = 1:NSIM
                    pointall(ii2) = pirf0{ii2}(il,ij,ik)/pirf0{ii2}(1,ik,ik);
                end
                crb(il,ij,ik) = prctile(pointall,PCTL(ii)*100);
            end
        end
    end
    CFBD{ii} = crb;    
end

IRF = pirf0{NSIM+1};
for ii = 1:nvar
    IRF(:,:,ii) = IRF(:,:,ii)/IRF(1,ii,ii);
end    

end
    