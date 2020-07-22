 % 
 % First cut at developing a ray integrator for EM and other waves:
 % (JFM May 13, May 20, May 21, 2020)
 %
 
 global cnst
 
 cnst = initCnst;           % will put more things in initCnst...
 
 % get the hydro profile
 dracoFile = "draco_EPsph_JFM.mat";
 tslice = 13;                % <- set the time slice here
 
 if ~exist('rayGd','var')
     disp("loading hydro...")
     rayGd = importDracoGrid(dracoFile,tslice);   
     disp("done loading hydro")
 end
 
 % update if the time slice has changed too
 if rayGd.iTime ~= tslice
     disp("updating hydro...")
     rayGd = importDracoGrid(dracoFile,tslice);
     disp("done updating hydro")
 end

 
 %
 % initialize and create a launch list for an EM wave ray bundle
 %

 launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
 launchList.mode = 'forward';        % Could be backward also (neg omega?).
 launchList.nrays = 20;
 launchList.frequency = cnst.omega0; % 1/sec
 launchList.focalPt = [-400,0];      % microns
 launchList.spot = struct('type','SG4','diameter',450);
 angle = 23.0; % (degrees) think about how we want to define the angle
 launchList.centroid = -[cosd(angle),sind(angle)]; % unit vector in
                                                   % direction of
                                                   % beam propagation
                                 % propagation
 launchList.translate = 5.0e3;   % distance in um from focus to
                                 % translate so that we are sure to
                                 % be far enough away to start

 % Create a ray bundle and bring the ICs into the sim volume. Then
 % we can start to integrate
 rayBundle = makeRayBundle(launchList);
 rayBundle = moveToDomain(rayBundle,rayGd);
 
 %
 % I'll probably move the integration below to a function
 %
 
 %
 % Integrate all the rays in the bundle one at a time
 %
 for rayIdx = 1:rayBundle.nrays; 
 
     % adjust intial condition to sit on the dispersion surface
     %
     omega = rayBundle.frequency;       % 1/s
     omega_ps = 1.e-12*omega;           % 1/ps
     kVac = 1.e-6*abs(omega)/cnst.c;    % 1/um
     x0 = rayBundle.rayICs(:,rayIdx)';  % row vector here
     localNe = 10.^interpOnTraj('valsNe',[x0 x0],rayGd);
     nc = rayBundle.nc;                 % 1/cm3
     k0Mag = kVac*sqrt(1-localNe/nc);   % 1/um
     kdir = rayBundle.direction;        % row vector
     k0 = k0Mag*kdir;                   % initial k (row) vector
                                        % (1/um)
     ray0 = [x0, k0]';                  % initial condition (column vector)
                                        % in phase space for ode integrator
 
     % The time integration should be done in short spurts so that we
     % can test to see if we need to do another one, or perhaps deal
     % with a caustic etc.
 
     tspan = [0 8];     % (ps) time range to solve over 
                        %   in future: tspan=[rays.time rays.time+tintv];  
 
     % Now we can integrate over the given time span and see what it
     % looks like
     %
     [tr,yr] = ode45(@(t,y) odeEmRayFun(t,y,omega_ps,rayGd),tspan, ...
                     ray0);
     
      % attach solution to rayBundle structure
     rayBundle.trajs{rayIdx} = [tr,yr];
 
     % might want to think about appending if this is not the first
     % integration (i.e., a restart)
 end
 
 
 % Let's plot it on the density so that we can visualize the ode
 % solution
 
 debug.density = true;
 
 if debug.density   % with the straight line
     figure(2)
     clf
     pcolor(rayGd.zSblock,rayGd.rSblock,rayGd.neblock);
     shading interp
     colorbar
     hold on
 
     for iplt = 1:rayBundle.nrays
         plot(rayBundle.trajs{iplt}(:,2),rayBundle.trajs{iplt}(:, ...
                                                           3),'r')
     end
     
     title('Log10 of electron density (1/cm3)');
     xlabel('z in um');
     ylabel('r in um');          

     %legend('str line');
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",rayGd.time))

     
     hold off
 end

 debug.temperature = true;
 
 if debug.temperature   % with the straight line
     figure(3)
     clf
     pcolor(rayGd.zSblock,rayGd.rSblock,rayGd.teblock);
     colormap('hot')
     shading interp
     colorbar
     hold on
 
     for iplt = 1:rayBundle.nrays
         plot(rayBundle.trajs{iplt}(:,2),rayBundle.trajs{iplt}(:, ...
                                                           3),'k')
     end

     % add nc
     qcrit = log10(rayBundle.nc);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qcrit ...
                         qcrit],'k:')     
     
     % add nc/4 too
     qcrit = log10(rayBundle.nc/4);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qcrit ...
                         qcrit],'k:')
     % add nc/10
     n10 = log10(rayBundle.nc/10);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[n10 ...
                         n10],'k:')
     
     title('electron temperature (eV)');
     xlabel('z in um');
     ylabel('r in um');     
     
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",rayGd.time))
     
     %legend('str line');
     hold off
 end
 
 
 %
 % Now let's integrate Raman light back from the detector
 %

 ramanList.type = 'RamanAtDetector';      % trigger for 'makeRayBundle'
 ramanList.mode = 'backward';   % integrate backwards (neg omega).
 ramanList.nrays = 18;

 % derive frequency from detected wavelength in nm
 %detectWavl = (700.0)*1.e-9;    % meters
 detectWavl = (650.0)*1.e-9;    % meters
 %detectWavl = (600.0)*1.e-9;    % meters
 
 ramanList.frequency = cnst.twopi*cnst.c/detectWavl;  % 1/sec

 ramanList.focalPt = [-400,0];      % microns
 ramanList.spot = struct('type','SG4','diameter',450);
 angle = 23.0; % (degrees) think about how we want to define the angle
 ramanList.centroid = [cosd(angle),sind(angle)];   % unit vector in
                                                   % direction of
                                                   % beam propagation
                                 % propagation
 ramanList.translate = 5.0e3;   % distance in um from focus to
                                 % translate so that we are sure to
                                 % be far enough away to start

 % Create a ray bundle and bring the ICs into the sim volume. Then
 % we can start to integrate
 ramanBundle = makeRayBundle(ramanList);
 ramanBundle = moveToDomain(ramanBundle,rayGd);


 % Integrate all the rays in the bundle one at a time
 %
 for rayIdx = 1:ramanBundle.nrays; 
 
 %rayIdx = 10;
     % adjust intial condition to sit on the dispersion surface
     %
     omega = ramanBundle.frequency;       % 1/s
     omega_ps = 1.e-12*omega;             % 1/ps
     kVac = 1.e-6*abs(omega)/cnst.c;      % 1/um
     x0 = ramanBundle.rayICs(:,rayIdx)';  % row vector here
     localNe = 10.^interpOnTraj('valsNe',[x0 x0],rayGd);
     nc = ramanBundle.nc;                 % 1/cm3
     k0Mag = kVac*sqrt(1-localNe/nc);     % 1/um
     kdir = ramanBundle.direction;        % unit row vector
     k0 = k0Mag*kdir;                     % initial k (row) vector
                                          % (1/um)
     ray0 = [x0, k0]';                    % initial condition (column vector)
                                          % in phase space for ode integrator
 
     % The time integration should be done in short spurts so that we
     % can test to see if we need to do another one, or perhaps deal
     % with a caustic etc.
 
     tspan = [0 7];     % (ps) time range to solve over 
                        %   in future: tspan=[rays.time rays.time+tintv];  
 
     % Now we can integrate over the given time span and see what it
     % looks like
     %
     [tr,yr] = ode45(@(t,y) odeEmRayFun(t,y,omega_ps,rayGd),tspan, ...
                     ray0);
     
      % attach solution to rayBundle structure
     ramanBundle.trajs{rayIdx} = [tr,yr];
 
     % might want to think about appending if this is not the first
     % integration (i.e., a restart)
 end
 
 debug.srsFromDetector = false;
 
 if debug.srsFromDetector
 
     figure(2)
     hold on

     for iplt = 1:ramanBundle.nrays
     
         plot(ramanBundle.trajs{iplt}(:,2),ramanBundle.trajs{iplt}(:, ...
                                                           3),'g-')
         % add a contour for the approx. Raman density
         omega_epw = rayBundle.frequency-abs(ramanBundle.frequency);
         ne_epw = (omega_epw/cnst.wpe)^2;
         cLine = log10(ne_epw);
         contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[cLine ...
                         cLine],'k')
         % add nc too
         qcrit = log10(rayBundle.nc);
         contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qcrit ...
                         qcrit],'k:')
         % add nc/4 too
         qcrit = log10(rayBundle.nc/4);
         contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qcrit ...
                         qcrit],'k:')
         % add nc/10
         n10 = log10(rayBundle.nc/10);
         contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[n10 ...
                         n10],'k:')     
     end
     
     hold off
 end
 
 
 %
 % Integrate Raman backscatter from an incident ray
 %
 %    the following should be moved into the makeRayBundle
 %    function...
 
 chosenRay = 18;
 ramanOnTraj = getRamanWavevectors(rayBundle.trajs{chosenRay}, ...
                                 rayBundle.frequency*1.e-12,rayGd);
 
 % Raman light
 yRaman = [ramanOnTraj{5} ramanOnTraj{6}];
 omRaman = ramanOnTraj{4};
 srsOnRayBundle.trajs = cell(1,numel(omRaman));
 
 % Raman EPW
 yEPW = [ramanOnTraj{2} ramanOnTraj{3}];
 omEPW = ramanOnTraj{1};
 epwOnRayBundle.trajs = cell(1,numel(omEPW)); 
 
 tspan = [0,2];
 tspanEPW = [0,6];

 
 % only one raman event 
 %idx = 13;
 
 % light
 for idx=1:5
     omega_ps = omRaman(idx);
     ray0 = yRaman(idx,:)';
     [trbs,yrbs] = ode45(@(t,y) odeEmRayFun(t,y,omega_ps,rayGd),tspan, ...
                     ray0);
     srsOnRayBundle.trajs{idx} = [trbs,yrbs];
 end

 % EPW
 for idx=1:5
     omega_ps = omEPW(idx);
     ray0 = yEPW(idx,:)';
     [tepw,yepw] = ode45(@(t,y) odeLwRayFun(t,y,omega_ps,rayGd),tspanEPW, ...
                     ray0);
     epwOnRayBundle.trajs{idx} = [tepw,yepw];
 end
 
 
 debug.srsOnTraj = true;
 
 if debug.srsOnTraj
   figure(2)
   hold on
 
   %   iplt = idx;
   for iplt = 1:5
       plot(epwOnRayBundle.trajs{iplt}(1,2), ...
            epwOnRayBundle.trajs{iplt}(1,3), 'oc')
       plot(epwOnRayBundle.trajs{iplt}(:,2), ...
            epwOnRayBundle.trajs{iplt}(:,3), '-c')
       %       plot(srsOnRayBundle.trajs{iplt}(:,2), ...
       %    srsOnRayBundle.trajs{iplt}(:,3), 'om')

   end
   
   for iplt = 1:5
       plot(srsOnRayBundle.trajs{iplt}(1,2), ...
            srsOnRayBundle.trajs{iplt}(1,3), 'om')
       plot(srsOnRayBundle.trajs{iplt}(:,2), ...
            srsOnRayBundle.trajs{iplt}(:,3), '-m')
       %       plot(srsOnRayBundle.trajs{iplt}(:,2), ...
       %    srsOnRayBundle.trajs{iplt}(:,3), 'om')

   end
 
   hold off
 end
 
 % 
 % Local functions follow:
 %

 function [t,y] = odeStrLine(tspan,y0,npts)
 % just evolve the ray trajectory in a straight line segment
    global cnst
    
    t = linspace(tspan(1),tspan(2),npts)'; % ps
    ttmp = [t t];
    k0 = y0(3:4)';
    kdir = k0/sqrt(dot(k0,k0)); 
    x0 =y0(1:2)';
    deltaX = cnst.c*(1.e-6)*kdir;        % [vz vr]
    x = repmat(x0,npts,1) + ttmp.*repmat(deltaX,npts,1);   % c
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
     
    global cnst

    clum  = (cnst.c)*(1.e-6); % speed of light in microns/ps
    ln10  = cnst.ln10;
    twopi = cnst.twopi;

    % need these
    lambdaum = twopi*clum/abs(omega_ps); % vac wavelength microns
    kVac = abs(omega_ps)/clum;           % vacuum wavenumber 1/um
    nc = 1.1e21/lambdaum^2;              % crit density in 1/cm^3

    x = y(1:2);     % current position at phase space point y
    kVec = y(3:4);  % current ray wavevector at phase space point y
    
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
    
 end

 
 function dydt = odeLwRayFun(t,y,omega_ps,rayGd)
 % RHS for our ray ode for Langmuir waves
 %  output dydt is a column vector:
 %  dydt = [dz/dt,dr/dt,dk_z/dt,dk_r/dt]'
 
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

