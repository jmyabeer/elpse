 % 
 % First cut at developing a ray integrator for EM waves:
 %   We need an example to see how ODE45 works
 %   particularly with our interpolating function
 % (JFM May 13, 2020)
 
 % To start you need to import "dracoGrid"
 
 % need to import our constants and plasma parameters; probably
 % want to have the grid in a global structure too
 
 global cnst rayGd miscPass
 
 cnst = initCnst;           % will put more things in initCnst
 
 lambdaum = 0.351;          % microns
 lambda0 = 0.351e-6;        % (m) probably this should live in the ray
                            % structure
 kVac = 2*cnst.pi/lambdaum; % inverse microns
 nc = 1.1e21/lambdaum^2;    % cm^-3 figure out where to put this
 omega0 = 2.*cnst.pi*cnst.c/lambda0; % (1/sec) this too
 
 % Initial condition
 %   given the plasma parameters at the position x and the
 %   frequency of the wave we need to solve the dispersion relation
 %   to obtain the magnitude of the wave vector. We can then pick
 %   the direction (and polarization) as we wish and integrate
 
 omegaRay = omega0;   % (1/sec) this should be held in a ray
                      % structure and extracted as needed
 
 % just to test...
 z0 = 500;            % (microns) the initial position in conf space will
 r0 = 1000;            % come from somewhere
 
 x0 = [z0, r0];       % conf space position in terms of inital z and r
 
 % Get the plasma density at this point
 %    Let's assume we have the block data read in from "dracoGrid.m"
 
 localNe = 10.^interpOnTraj('valsNe',[x0 x0]);   % ne on trajectory
                                                 % duplicate since
                                                 % we don't know k0
 
 k0Mag = kVac*sqrt(1-localNe/nc);          
                      % At some point do this in more generality
                      % assuming a given dispersion matrix. See for
                      % example T =
                      % dispertok(0,y0(:,i),y0(:,i),0,'Msw')
 
 kdir = [-1,-1];      % initial direction for the ray [zcmpt,rcmpt]
 kdir = kdir/sqrt(dot(kdir,kdir));   % normalized
 k0 = k0Mag*kdir;     % initial k vector (1/um)
 
 ray0 = [x0, k0]';    % initial condition (column vector) in phase space
 
 % The time integration should be done in short spurts so that we
 % can test to see if we need to do another one, or perhaps deal
 % with a caustic etc.
 
 tspan = [0 6.0];   % (ps) time range to solve over 
                    %     in future: tspan=[rays.time rays.time+tintv];  
 
 % Now we can integrate over the given time span and see what it
 % looks like
 
 % Let's just do a straight line segment for now and plot it on the
 % density so that we can compare it with the ode solution
  
 [tl,yl] = odeStrLine(tspan,ray0,50);
 
 rayNe = 10.^interpOnTraj('valsNe',yl);   % ne on trajectory
 rayTe = interpOnTraj('valsTe',yl);
 rayTi = interpOnTraj('valsTi',yl);
 % and the derivatives...
 rayDLogNedz = interpOnTraj('valsDLogNedz',yl);
 rayDLogNedr = interpOnTraj('valsDLogNedr',yl);
 
 % Let's write a function to plot hydro quantities defined on a
 % trajectory (later)
 
 debug.derivs = true;
 
 if debug.derivs
     figure(1)
     clf
     plot(tl,rayDLogNedz)
     hold on
     plot(tl,rayDLogNedr)
     hold off
 end
 
 % NOW: go back and do the same trajectory with ode45
 
 [tr,yr] = ode45(@odeEmRayFun,tspan,ray0);
 
 
 % Integrate a LW trajectory
 
 % change this into a look over the rayBundle after it has been
 % moved into the simulation domain
 
 % move each ray into the domain and adjust the rayICs
 % appropriately
 
 % then integrate each ray for a given amount of time
 
 z0 = 400;            % (microns) the initial position in conf space will
 r0 = 300;             % come from somewhere
 
 x0 = [z0, r0];        % conf space position in terms of inital z and r
 
 localNe = 10.^interpOnTraj('valsNe',[x0 x0]);   % ne on trajectory
 localTe = interpOnTraj('valsTe',[x0 x0]);       % Te on trajectory

 lamDebye = 1.e4*(cnst.lamDebye)*sqrt(localTe/localNe); % microns
 wpe = cnst.wpe*sqrt(localNe)*1.e-12;   % rad/ps

 % looks like we have an error in kLw's magnitude
 kLw = 0.2/lamDebye;   % initial wave vector magnitude
 omegaLw = wpe*sqrt(1+3.0*(0.2^2));
 miscPass.omega = omegaLw;

 kdir = [-1,0];      % initial direction for the ray [zcmpt,rcmpt]
 kdir = kdir/sqrt(dot(kdir,kdir));   % normalized
 k0 = kLw*kdir;     % initial k vector (1/um)

 
 ray0Lw = [x0, k0]';   % initial condition (column vector) in phase space
 tspanLw = [0 60.0];    % (ps) time range to solve over 


 % now integrate the trajectory
 
 [tLw,yLw] = ode45(@odeLwRayFun,tspanLw,ray0Lw);
 
 %             rhs,  tspan, y0 (column vector)
 % each row of y is solution for a given time
 % t is a column vector of times
 
 % Plot to see what it looks like:
 %   The plotting stuff should be done using a seperate function
 %   (later)
 
 % the following needs to be reworked
 
 debug.density = true;
 
 if debug.density   % with the straight line
     figure(2)
     clf
     pcolor(rayGd.zSblock,rayGd.rSblock,rayGd.neblock);
     shading interp
     colorbar
     hold on
 
     plot(yl(:,1),yl(:,2),'k-')
     plot(yr(:,1),yr(:,2),'r')
     
     plot(yLw(:,1),yLw(:,2),'g-')

     title('Solution of ray equations with ODE45');
     xlabel('z in um');
     ylabel('r in um');          
     %legend('str line');
     hold off
 end

 
% Local functions:
 
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

 
 function varOnTraj = interpOnTraj(varToGet,y)
 % Interpolates hydro variables onto a ray segment
    global rayGd

    x = y(:,1:2);
    [ti,bc] = pointLocation(rayGd.DT,x);     % get Delauney
                                             % triangles
    switch varToGet
      case 'valsNe'
        % to do - check to see if the variable exists
        triVals = rayGd.valsNe(rayGd.DT(ti,:));
      case 'valsDLogNedz'
        triVals = rayGd.valsDLogNedz(rayGd.DT(ti,:));
      case 'valsDLogNedr'
        triVals = rayGd.valsDLogNedr(rayGd.DT(ti,:));
      case 'valsTe'
        triVals = rayGd.valsTe(rayGd.DT(ti,:));
      case 'valsTi'
        triVals = rayGd.valsTi(rayGd.DT(ti,:));
      otherwise
        error("No matching hydro variable found in interpOnTraj");
    end
    
    varOnTraj = dot(bc',triVals')';
 end

 
 function dydt = odeEmRayFun(t,y)
 % RHS for our ray ode for electromagnetic waves
 %  output dydt is a column vector:
 %  dydt = [dz/dt,dr/dt,dk_z/dt,dk_r/dt]'
 
    global cnst rayGd

    clum  = (cnst.c)*(1.e-6); % speed of light in microns/ps
    ln10  = cnst.ln10;
    twopi = cnst.twopi;

    % these need to be stored somewhere better:    
    lambdaum = 0.351;               % vac wavelength microns
    kVac = 2*(cnst.pi)/lambdaum;    % vacuum wavenumber 1/um
    omega = twopi*clum/lambdaum;    % time in picoseconds
    nc = 1.1e21/lambdaum^2;         % crit density 1/cm^3

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
    
    dzdt = clum*kVec(1)/kVac;
    drdt = clum*kVec(2)/kVac;
    
    dkzdt = 0;  % for now  (free space propagation)
    dkrdt = 0;  % for now

    % overwrite 
    dkzdt = -0.5*ln10*omega*netonc*dLogNedz;
    dkrdt = -0.5*ln10*omega*netonc*dLogNedr;
    
    dydt = [dzdt,drdt,dkzdt,dkrdt]'; % column vector 
    
 end

 
 function dydt = odeLwRayFun(t,y)
 % RHS for our ray ode for Langmuir waves
 %  output dydt is a column vector:
 %  dydt = [dz/dt,dr/dt,dk_z/dt,dk_r/dt]'
 
    global cnst rayGd miscPass

    clum  = (cnst.c)*(1.e-6); % speed of light in microns/ps
    ln10  = cnst.ln10;
    twopi = cnst.twopi;

    % these need to be stored somewhere better:    
    lambdaum = 0.351;               % vac wavelength microns
    ncLaser = 1.1e21/lambdaum^2;              % crit density 1/cm^3
    
    kVac = 2*(cnst.pi)/lambdaum;    % vacuum wavenumber 1/um
    laserOmega = twopi*clum/lambdaum;    % time in picoseconds

    omega = miscPass.omega;

    x = y(1:2);     % current position at phase space point y
    kVec = y(3:4);  % current ray wavevector at phase space point y
    k2 = dot(kVec,kVec);
    
    % interpolation for current position
    [ti,bc] = pointLocation(rayGd.DT,x');  % Delauney triangles    
    
    % density
    triValNe = rayGd.valsNe(rayGd.DT(ti,:));
    logNe = dot(bc',triValNe')';           % log10 of electron density
%    disp(logNe)  % debugging
    netoncL = 10^(logNe)/ncLaser;    
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
    
    dkzdt = 0;  % for now  (free space propagation)
    dkrdt = 0;  % for now

    % need to debug
    
    % overwrite 
    dkzdt = -wpe^2/(2.0*omega)*ln10*dLogNedz - (3/2)*(k2*vTe2/omega)*dLnTedz;
    dkrdt = -wpe^2/(2.0*omega)*ln10*dLogNedr - (3/2)*(k2*vTe2/omega)*dLnTedr;
    
    dydt = [dzdt,drdt,dkzdt,dkrdt]'; % column vector 
    
 end

