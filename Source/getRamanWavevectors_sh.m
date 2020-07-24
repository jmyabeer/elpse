function [srsBundle,epwBundle] = getRamanWavevectors_sh(traj,omega_ps, ...
                                                       rayGd, angle, ...
                                                       itmax, ...
                                                       landauC,inboundOnly)    
% function [srsBundle,epwBundle] =
% getRamanWavectors_sh(traj,omega_ps,rayGd,angle,itmax,landauC,inboundOnly)
%
%         traj - the trajectory of a single EM ray
%                = time z r kz kr, where each is a column vector
%     omega_ps - the (scalar) frequency of the ray trajectory (rads/ps)
%        rayGd - output struct from reading hydrodynamics
%        angle - angle of SRS to be made in degrees
%        itmax - unused here?
%      landauC - value to use for the Landau cutoff - e.g. 0.3
%  inboundOnly - do not make Raman after the ray has encountered nc/4  


    global cnst   
    cumps = (cnst.c)*1.e-6;       % in um/ps

    % need to know the frequency of the ray
    kvac = omega_ps/cumps;        % (omega_ps is a scalar) inverse microns
    lambdaVac = cnst.twopi/kvac;  % microns
    nc = 1.1e21/lambdaVac^2;      % cm^-3    
    
    [nrows ncols] = size(traj);
    
    % treat omega_ps as a scalar for now
    if isscalar(omega_ps)
        omega = repmat(omega_ps,nrows,1);
    else
        error('omega_ps must be a scalar')
        % at some point we'll allow it to vary over a trajectory
    end
    
    % default number of iterations
    if ~exist('itmax','var')
        itmax = 5;
    end    

    % default Landau cutoff
    if ~exist('landauC','var')
        landauC = 0.3;   % maximum allowable k*lambdaDebye 
    end    

    % default inboundOnly is true (no SRS after nc/4 is reached)
    if ~exist('inboundOnly','var')
        inboundOnly = true;  
    end    
    
    % For each point on ray we want a density and temperature    
    neList = interpOnTraj('valsNe',traj,rayGd);
    teList = interpOnTraj('valsTe',traj,rayGd);   % eV

    % let's first assume backscatter and then generalize it later
    if ncols == 4
        xLoc  = traj(:,1:2);   % (Z, R)    (nrowsx2)
        kvecs = traj(:,3:4);   % (kZ,kR)   (nrowsx2)
    elseif ncols == 5
        tt = traj(:,1);        % times
        xLoc  = traj(:,2:3);   % (Z, R)
        kvecs = traj(:,4:5);   % (kZ,kR)
    end

    % incident light wavevector magnitude for each point
    k0 = sqrt(kvecs(:,1).^2 + kvecs(:,2).^2);  % um^-1 (nrowsx1)
    
    k0Hat = kvecs./[k0 k0];   % unit vector in the direction of k0     
    
    % plasma paramters on the ray points
    ne = 10.^(neList);     % cm^-3
    neTonc = ne/nc;
    wpeTow0 = sqrt(neTonc);   
    lamDeb2 = (cnst.lamDebye)^2*teList./ne;    % cm^2
    lamDebum2 = 1.e8*lamDeb2;                  % um^2
    
    %
    % We might want to exclude Raman on the outbound trajectory,
    % i.e. it is considered to be outbound if the ray has already
    % encountered quarter critical density

    % let's see if this is the case
    %
    outBound = find(neTonc >= 0.25); % indices for outbound points

    if isempty(outBound)
        lastInboundIdx = numel(neTonc);
    else
        lastInboundIdx = outBound(1);
    end
    
    % only want a Raman event if the solution for k is real - we
    % can remove those with imaginary parts when we are
    % done. i.e. it does no harm to compute them so long as we
    % REMEMBER this
    
    % Raman EPW wavevector (k) magnitude (backscatter); iterate klamD
    % corrections
    %
    
    rotMat = [cosd(angle),-sind(angle);sind(angle),cosd(angle)];
    
    % kSRSHat is rotated counterclockwise by "angle" degrees from k0Hat
    for i = 1:nrows 
        columnVec = rotMat*[k0Hat(i,1);k0Hat(i,2)];
        kSRSHat(i,1) = columnVec(1);
        kSRSHat(i,2) = columnVec(2);
    end
    
    kSRSMag = k0.*sqrt(1-2.*sqrt(neTonc));                 
    kSRS = kSRSMag.*kSRSHat;
    kLangmuir = kvecs-kSRS;
    k = sqrt(kLangmuir(:,1).^2 + kLangmuir(:,2).^2);
    k2lam2 = (k.^2).*lamDebum2;
    
    % can't have Raman just anywhere ...
    %
    % - This test wasn't good enough. We should test on
    %   kSRSMag instead of k (JFM)
    %
    
    isRealK = imag(kSRSMag) == 0;               % logical vectors
    isBelowLandau = abs(k2lam2) <= landauC^2;
    isInbound = (1:numel(k))' <= lastInboundIdx;   % think how to do this

    if inboundOnly
        goodidxs = find(isRealK & isBelowLandau & isInbound);
    else
        goodidxs = find(isRealK & isBelowLandau);
    end    
    
    % Get the times
    %
    if exist('tt','var')
        ttgood = tt(goodidxs);
    else
        ttgood = zeros(1,numel(goodidxs));
    end
    
    % probably want to change to below to be computed just on the
    %    good indices
    freqEPW = omega_ps.*wpeTow0.*sqrt(1+3*k2lam2);        % ps^-1
    %freqEPW = omega_ps.*wpeTow0.*(1+3/2*k2lam2);        % ps^-1    
    freqSRS2 = (omega_ps.*wpeTow0).^2 + cumps^2*kSRS.^2;  % ps^-2
    freqSRS = sqrt(freqSRS2);                             % ps^-1
    vacWavlSRS = cnst.twopi*cumps./freqSRS;  % microns
        
    %
    % Langmuir wave first
    %
    xs = xLoc(goodidxs,:);         % row vector (z r)
    ks = kLangmuir(goodidxs,:);    % row vector (kz kr)
    
    epwBundle.name = 'plasmaWave';
    epwBundle.type = 'EPW';
    epwBundle.mode = 'forward';
    epwBundle.nrays = numel(goodidxs);
    epwBundle.halt = zeros(1,numel(goodidxs));
    epwBundle.frequency = 1.e12*freqEPW(goodidxs)';   % (row vec) s^-1
    epwBundle.trajs = cell(1,numel(goodidxs));
    
    % put in an intial condition for each trajectory
    for idx = 1:numel(goodidxs)
        % needs to be a row vector: t, z, r, kz, kr
        epwBundle.trajs{idx} = [ttgood(idx) xs(idx,:) ks(idx,:)];
    end
    
    %
    % then Raman light
    %
    xs = xLoc(goodidxs,:);         % row vector (z r)
    ks = kSRS(goodidxs,:);    % row vector (kz kr)
    
    srsBundle.name = 'srsLight';
    srsBundle.type = 'EM';
    srsBundle.mode = 'forward';
    srsBundle.nrays = numel(goodidxs);
    srsBundle.halt = zeros(1,numel(goodidxs));
    srsBundle.frequency = 1.e12*freqSRS(goodidxs)'; % row vec 1/sec
    lam0 = (cnst.twopi)*(cnst.cumps)./freqSRS(goodidxs)';
    srsBundle.nc = 1.1e21./lam0.^2; % cm^-3 for use in pushBundle()
    srsBundle.trajs = cell(1,numel(goodidxs));
    
    % put in an intial condition for each trajectory
    for idx = 1:numel(goodidxs)
        % needs to be a row vector: t, z, r, kz, kr
        srsBundle.trajs{idx} = [ttgood(idx) xs(idx,:) ks(idx,:)];
    end

end

