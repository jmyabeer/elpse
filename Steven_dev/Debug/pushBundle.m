
function bundleOut = pushBundle(rayBundle,rayGd,tStep,margin,npts)
% 
% First cut at developing a ray integrator for EM and other waves:
% (JFM May 28, 2020)
%
% Perhaps trajIsNew should be replaced with pushFromBoundary? This
%  might help with the problems associated with integrating 
%  trajectories that could potentially leave the domain, and
%  where we currently have problems with "margin" (JFM 21 July, 2020).
%
 
 global cnst
 
 if ~exist('cnst','var')
     cnst = initCnst;              % will put more things in initCnst...
 end

 if ~exist('margin','var')
     margin = [150 150 150 150];   % margin around domain in microns
 end

 if ~exist('npts','var')
     npts = 2000;                  % points to use in str line check
 end
 
 %
 % Integrate all the rays in the bundle one at a time
 %
 
 for rayIdx = 1:rayBundle.nrays; 
 
     % check to see if we need to move this ray
     if rayBundle.halt(rayIdx)
         continue  % skip to next ray (iteration of loop)
     end
     
     trajIsNew = isempty(rayBundle.trajs{rayIdx}); 
     waveType = rayBundle.type;
     
     if trajIsNew    % need IC's for trajectory
         if waveType == 'EM'
             x0 = rayBundle.rayICs(:,rayIdx)';  % position x0 is row vector here
             % modify arguments to "toDispSurface()" to only what is
             % needed
             %
             k0 = toDispSurface(x0,rayBundle,rayIdx,rayGd);         
             tPrev = 0;                         % start from t=0
         else
             error('Currently new trajectories are for EM waves only')
         end         
     else
         lastInfo = rayBundle.trajs{rayIdx}(end,:); % use the last time info
         x0 = lastInfo(2:3);   % position row vector
         k0 = lastInfo(4:5);   % wave vector row vector
         tPrev = lastInfo(1);
     end

     % TO DO: deal with the case where rayBundle.frequency is a row vector
     %
     omega_ps = 1.e-12*rayBundle.frequency(rayIdx);
     dirn = sign(omega_ps);             % see which direction we're going
     tStop = tPrev+tStep;
     tSpan = [tPrev, tStop];     
     ray0 = [x0, k0]';                  % initial condition (column vector)
                                        % in phase space for ode integrator

% $$$      % Check to see if we could fall outside of the interpolation
% $$$      % domain...
% $$$      % need to pass the wave type to odeStrLine (i.e., for EPW)
% $$$      %
% $$$      [tTest,yTest] = odeStrLine(tSpan,ray0,npts,dirn,waveType,rayGd);
% $$$ 
% $$$      if trajIsNew
% $$$          indom = inDomain(yTest(:,1:2),rayGd);  % column vec of
% $$$                                                 % logicals
% $$$      else
% $$$          indom = inDomain(yTest(:,1:2),rayGd,[10 10 10 10]);  
% $$$      end
% $$$      
% $$$      badIdxs = find(~indom);
% $$$      
% $$$      % update upper time limit if necessary
% $$$      if any(badIdxs)
% $$$          dts = diff(tSpan)/npts;
% $$$          badTimes = tTest(badIdxs);
% $$$          tStop = max(tPrev,badTimes(1)-dts);  % happens if all
% $$$                                               % times are bad 
% $$$          tSpan = [tPrev,tStop]; % adjust upper limit
% $$$      end
     
     % Provided we have a good integration time we can integrate
     % over the given time span and see what it looks like
     %
     if tStop > tPrev
         switch waveType
           case 'EM'
             [tr,yr] = ode45(@(t,y) odeEmRayFun(t,y,omega_ps, ...
                                                rayGd),tSpan,ray0);
           case 'EPW'
             [tr,yr] = ode45(@(t,y) odeLwRayFun(t,y,omega_ps, ...
                                                rayGd),tSpan,ray0);           
           otherwise
             exit('invalid waveType')
         end
         
         newTraj = [tr,yr];   % because otherwise ode45 will return
                              % a struct         
         % need to check and set halt flags if we are too close to
         % theboundary at the last time
         %
         withinMargin = inDomain(yr(end,1:2),rayGd,margin);
         if any(~withinMargin)
             rayBundle.halt(rayIdx) = 1;   % set halt flag for
                                           % out-of-margin
         end
         % attach solution to rayBundle structure
         %  have two cases here... 
         if trajIsNew
             rayBundle.trajs{rayIdx} = newTraj;
         else
             % append to existing
             oldTraj = rayBundle.trajs{rayIdx}(1:end-1,:); % drop IC
             rayBundle.trajs{rayIdx} = [oldTraj; newTraj];
             %             rayBundle.trajs{rayIdx} = newTraj;   % debug
         end
       end  % the if tSpan..
          
 end   % for loop over rays
 
 bundleOut = rayBundle; 
 
end 
 
 % 
 % Local functions follow:
 %

 function k0 = toDispSurface(x0,rayBundle,rayIdx,rayGd)
 % TODISPSURFACE - to do later
 %   
 % TO DO: need this function to work on EPW also
     
     global cnst
     
     % get magnitude of local wavevector
     omega = rayBundle.frequency(rayIdx);       % 1/s
     % TO DO: need to deal with the case that rayBundle.frequency
     % is a row vector
     %
     omega_ps = 1.e-12*omega;           % 1/ps
     % we assume light - need to fix to be more general
     %    we can use rayBundle.type for this...
     kVac = 1.e-6*abs(omega)/cnst.c;    % 1/um

     localNe = 10.^interpOnTraj('valsNe',[x0 x0],rayGd);
     nc = rayBundle.nc(rayIdx);         % 1/cm3
     k0Mag = kVac*sqrt(1-localNe/nc);   % local wavenumber in
                                        % 1/um
     kdir = rayBundle.direction;        % initial direction
                                        % (unit row vector)
     k0 = k0Mag*kdir;                   % initial k (row) vector
                                        % (1/um)
 end 
 
 
 function inRange = inDomain(x,rayGd,margin)
 % INDOMAIN - given a column vector of positions x, return a list
 % of logicals showing if they are in the domain (given by rayGd)
 % or not
 %   
     
    if ~exist('margin','var')
        margin = [0 0 0 0];    % margin in zmin zmaz rmin rmax (um)
    end
    
    rayZ = x(:,1);
    rayR = x(:,2);
    
    zmin = rayGd.domain(1) + margin(1);
    zmax = rayGd.domain(2) - margin(2);
    rmin = rayGd.domain(3) + margin(3);
    rmax = rayGd.domain(4) - margin(4);
    
    % could write as zmin = (rayGd.domain+margin)(1)
    
    inRrange = (rayR < rmax) & (rayR > rmin);
    inZrange = (rayZ < zmax) & (rayZ > zmin);
    inRange = inRrange & inZrange;
 end
 
 
 function [t,y] = odeStrLine(tspan,y0,npts,dir,waveType,rayGd)
 % just evolve the ray trajectory in a straight line segment at the
 % appropriate speed with npts and in the direction given by dir (= +1
 %    or -1) which accommodates the case when the ray is integrated
 %    backwards in time.
 %
 %    y0 is a column vector: [z r kz kr]'
 
    global cnst
    
    
    t = linspace(0,tspan(2)-tspan(1),npts)'; % time increment in ps
    ttmp = [t t];
    
    k0 = y0(3:4)';                 % turn into a row vector
    kdir = k0/sqrt(dot(k0,k0)); 
    x0 =y0(1:2)';                  % position row vector

    
    switch waveType
      case 'EM'
        speed = cnst.cumps;  % vacuum speed of light in um/ps
      case 'EPW'
        % interpOnTraj wants y0 to be a row vector
        %
        localTe = interpOnTraj('valsTe',y0',rayGd);   % eV
        speed = (cnst.vTe1eV)*sqrt(localTe);          % estimate using local Te      
      otherwise
        error('bad wave type in odeStrLine')
    end    
    
    deltaX = dir*speed*kdir;        % [vz vr]
    % x = x0 + t*v
    x = repmat(x0,npts,1) + ttmp.*repmat(deltaX,npts,1);   % speed
                                                           % needs
                                                           % to be
                                                           % in
                                                           % um/ps
    k = repmat(k0,npts,1);
    y = [x k];
 end

 
 function dydt = odeEmRayFun(t,y,omega_ps,rayGd)
 % RHS for our ray ode for electromagnetic waves
 %
 %     t        - time in ps
 %     y        - phase space point [z r kz kr] ('?)
 %     omega_ps - EM wave frequency in 1/ps
 %     rayGd    - grid struct 
 %  output dydt is a column vector:
 %  dydt = [dz/dt,dr/dt,dk_z/dt,dk_r/dt]'
 %

 %
 % Added a check to ensure that the position y(1:2) where
 % the RHS is to be evaluated is inside the domain where
 % the grid data is defined. Attempting to evaluate outside
 % of the domain will give "triangulation" errors, i.e.,
 % the pointLocation() function will give an error.
 %
 % If the point is outside the domain the function will
 % return zero for dydt in the second and third positions
 %
 % The same needs to be done for odeLwRayFun but I haven't
 % implemented it yet.
 %
 % (JFM 21 July, 2020)
     
    global cnst

    clum  = (cnst.c)*(1.e-6); % speed of light in microns/ps
    ln10  = cnst.ln10;
    twopi = cnst.twopi;

    % need these
    lambdaum = twopi*clum/abs(omega_ps); % vac wavelength microns
    kVac = abs(omega_ps)/clum;           % vacuum wavenumber 1/um
    nc = 1.1e21/lambdaum^2;              % crit density in 1/cm^3

    x = y(1:2);       % current position at phase space point y
    kVec = y(3:4);    % current ray wavevector at phase space point y

    goodPt = inDomain(x',rayGd);       
    
    if goodPt
        % interpolation for current position
        [ti,bc] = pointLocation(rayGd.DT,x');  % Delauney triangles    
        
        triValNe = rayGd.valsNe(rayGd.DT(ti,:));
        logNe = dot(bc',triValNe')';     % log10 of electron density
                                         %    disp(logNe)  % debugging
        netonc = 10^(logNe)/nc;
        
        triValDLogNedz = rayGd.valsDLogNedz(rayGd.DT(ti,:));
        dLogNedz = dot(bc',triValDLogNedz')';   % at phase space point
                                                %    disp(dLogNedz)

        triValDLogNedr = rayGd.valsDLogNedr(rayGd.DT(ti,:));
        dLogNedr = dot(bc',triValDLogNedr')';   % at phase space point
                                                %    disp(dLogNedr)
        
        dzdt = sign(omega_ps)*clum*kVec(1)/kVac;
        drdt = sign(omega_ps)*clum*kVec(2)/kVac;    
        dkzdt = -0.5*ln10*omega_ps*netonc*dLogNedz;
        dkrdt = -0.5*ln10*omega_ps*netonc*dLogNedr;
        
        dydt = [dzdt,drdt,dkzdt,dkrdt]'; % column vector 
    else
        % move ray in a straight line
        dzdt = sign(omega_ps)*clum*kVec(1)/kVac;
        drdt = sign(omega_ps)*clum*kVec(2)/kVac;    
        dkzdt = 0.0;
        dkrdt = 0.0;
        
        dydt = [dzdt,drdt,dkzdt,dkrdt]'; % column vector         
    end    
    
end

 
 function dydt = odeLwRayFun(t,y,omega_ps,rayGd)
 % RHS for our ray ode for Langmuir waves
 %  output dydt is a column vector:
 %  dydt = [dz/dt,dr/dt,dk_z/dt,dk_r/dt]'
     
 % NOTE: here omega_ps is a scalar (since we deal with a single ray)
 
    global cnst

    clum  = (cnst.c)*(1.e-6); % speed of light in microns/ps
    ln10  = cnst.ln10;
    twopi = cnst.twopi;

    omega = omega_ps;   % wave frequency in rads/ps

    x = y(1:2);     % current position at phase space point y
    kVec = y(3:4);  % current ray wavevector at phase space point y
    k2 = dot(kVec,kVec);
    
    % interpolation for current position
    [ti,bc] = pointLocation(rayGd.DT,x');  % Delauney triangles    
    
    % density
    triValNe = rayGd.valsNe(rayGd.DT(ti,:));
    logNe = dot(bc',triValNe')';           % log10 of electron density
%    disp(logNe)  % debugging
    wpe = (cnst.wpe)*sqrt(10^logNe)*1.e-12; % electron plasma
                                            % frequency (rad/ps)
%    disp(wpe)

    triValDLogNedz = rayGd.valsDLogNedz(rayGd.DT(ti,:));
    dLogNedz = dot(bc',triValDLogNedz')';   % at phase space point
%    disp(dLogNedz)

    triValDLogNedr = rayGd.valsDLogNedr(rayGd.DT(ti,:));
    dLogNedr = dot(bc',triValDLogNedr')';   % at phase space point
%    disp(dLogNedr)
    
    % derivatives of electron temperature
    triValTe = rayGd.valsTe(rayGd.DT(ti,:));
    Te = dot(bc',triValTe')';        % electron temperature in eV
    vTe2 =(cnst.vTe1eV*sqrt(Te))^2;  % square of electron thermal velocity (um/ps)^2
    
    triValDLnTedz = rayGd.valsDLnTedz(rayGd.DT(ti,:));
    dLnTedz = dot(bc',triValDLnTedz')';   % at phase space point
%    disp(dLnTedz)

    triValDLnTedr = rayGd.valsDLnTedr(rayGd.DT(ti,:));
    dLnTedr = dot(bc',triValDLnTedr')';   % at phase space point
%    disp(dLnTedr)
    
    dzdt = (3.0*vTe2/omega)*kVec(1);
    drdt = (3.0*vTe2/omega)*kVec(2);
    
    dkzdt = -wpe^2/(2.0*omega)*ln10*dLogNedz - (3/2)*(k2*vTe2/ ...
                                                      omega)*dLnTedz;
    
    dkrdt = -wpe^2/(2.0*omega)*ln10*dLogNedr - (3/2)*(k2*vTe2/ ...
                                                      omega)*dLnTedr;
    
    dydt = [dzdt,drdt,dkzdt,dkrdt]'; % column vector 
    
 end

