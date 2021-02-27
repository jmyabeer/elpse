 % 
 % Debugging the computation of SRS gain along a ray
 %
 %   Edited: JFM 21/OCT/2020
 
 
 path(path,'./Plotting')
 path(path,'./Source')
 %path(path,'./Steven_dev')
 
 global cnst
 
 % some initialization
 cnst = initCnst;             % will put more things in initCnst...
 pltIncl = initDefaultPlots;  % default plots to make
 pltIncl.temperature = true;  %   modify to add the elecron temp
 
 
 % get the hydro profile
 %
 dracoFile = "draco_EPsph_JFM.mat";

 % set the time slice here
 %
 tslice = 11;                

 % variables to import/define set here:
 %
 addVarFlag.ne = true;
 addVarFlag.dLogNedz = true;
 addVarFlag.dLogNedr = true;
 addVarFlag.Dmn = false;
 addVarFlag.te = true;
 addVarFlag.ti = false;
 addVarFlag.dLnTedz = true;
 addVarFlag.dLnTedr = true;
 addVarFlag.Vz = true;
 addVarFlag.Vr = true;
 addVarFlag.Zbar = true;
 addVarFlag.Zsqr = true;
 
 % TO DO: add a function called "importAnalyticGrid()" that allows
 % you to define analytic hydro profiles, e.g. linearly varying
 % density, constant temperature and so on. Everything will still
 % be stored in the rayGd structure though.
 %
 
 if ~exist('rayGd','var')
     disp("loading hydro...")
     rayGd = importDracoGrid(dracoFile,tslice,addVarFlag);   
     disp("done loading hydro")
 else
     disp("using exisiting hydro")
 end
 
 % update if the time slice has changed too
 if rayGd.iTime ~= tslice
     disp("updating hydro...")
     rayGd = importDracoGrid(dracoFile,tslice,addVarFlag);
     disp("done updating hydro")
 end

 
 %
 % initialize and create a launch list and ray bundle for an EM
 % wave - if one doesn't already exist
 %
 
 % At some point generalize makeRayBundle so that it can launch
 % EPWs (JFM)
 

 % Beam 1
 %
 if ~exist('rayBundleB1','var')
     
     launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
     launchList.mode = 'forward';        % Could be backward also (neg omega?).
     launchList.nrays = 20;
     launchList.frequency = cnst.omega0; % 1/sec
     % center of spherical target
     launchList.focalPt = [-400,0];      % microns
     launchList.spot = struct('type','SG8','diameter',700); 
     angle = 180+(-23.3); % (degres) is measured from "target norm"
     launchList.centroid = [cosd(angle),sind(angle)]; % unit vector in
                                                       % direction of
                                                       % beam propagation
                                                       % propagation
     launchList.translate = 5.0e3;   % distance in um from focus to
                                     % translate so that we are sure to
                                     % be far enough away to start
     % Create a ray bundle
     rayBundleB1 = makeRayBundle(launchList,rayGd);
     
     % give it a useful name
     rayBundleB1.name = 'Omega EP beam #1';
     
 end

 %
 %  Integrate bundles
 %

 % incident laser light - Beam 1
 
 % First push
 %
 tPush = 3.7; % ps
 
 rayBundleB1 = pushBundle(rayBundleB1,rayGd,tPush,[100 100 100 100]);
 
 % refine the ray push
 %
 nits = 35;
 
 for i=1:nits
     rayBundleB1 = pushBundle(rayBundleB1,rayGd,0.2);
 end

 % Halt any further integration of these rays
 %
 rayBundleB1.halt = setHaltAll(rayBundleB1);

 
 %
 %  Make some plots
 %
 
 % Let's plot it on the density so that we can visualize the ode
 % solutions
 
 makePlotList(pltIncl,rayGd);

 % add the incident beam and SRS 
 %
 figure(1)   % density plot
 hold on
 
 % add selected density contours (freq determined by beam #1)
 %
 addContourPlt(rayBundleB1,rayGd,'nc');
 addContourPlt(rayBundleB1,rayGd,'nc4');
 addContourPlt(rayBundleB1,rayGd,'nc10');
 
 % Add rays for incident beam #1
 %
 addBundlePlt(rayBundleB1,'k');

 
 %-------
 % Now let's try to compute the intensity
 %-------

 
 % choose a ray from the bundle
 testIdx = 3;
 
 % Manually put in the initial intensity for now:
 %   This should eventually be done by "makeRayBundle" using
 %   the laser spot shape as required by "makeLaunchList". i.e.,
 %   this information will be held in the ray bundle struct.
 %
 I0 = 5.e14;                                  % W/cm^2
 rayBundleB1.I0 = zeros(rayBundleB1.nrays,1); % one for each ray
 rayBundleB1.I0(testIdx) = I0;
 
 testTraj = rayBundleB1.trajs{testIdx};
 ncForRay = rayBundleB1.nc(testIdx);          % cm^-3
 
 gammaEM = interpOnTraj('gammaEM',testTraj,rayGd,ncForRay);
 
 % While we're at it, we can compute the path length variable too
 %    this is also stored in the ray bundle struct
 %
 rayBundleB1.path = computePathLength(rayBundleB1);
 
 % If desired, we could now plot the temporal absorption rate
 % (gammaEM) as a function of either time, or path length
 
 % As a function of time:
 %
 time = testTraj(:,1);
 tSamp = time;         % the times where gammaEM is known (sampled)
 
 
 % Can now solve the ode for intensity as a function of time along
 % the ray (ignoring ray divergence for now)
 %    We'll solve for the natural log of the intensity
 %
 
 % intial value
 %
 lnI0 = log(rayBundleB1.I0(testIdx));    % since we've decided to
                                         % solve for the log of I
 tspan = [time(1) time(end)];            % interval to solve the
                                         % ode over
 
 sol = ode45(@(t,y) dlnIdt(t,y,tSamp,gammaEM),tspan,lnI0);
 
 % now evaluate on our discrete trajectory
 %
 logI = deval(sol,time);
 
 % ... and assign the intensity to our ray struct.
 %
 rayBundleB1.I = cell(1,20);
 rayBundleB1.I{testIdx} = exp(logI);
 
 % you can now plot the intensity as a function of either time or
 % path length
 
 
 % Now we need to compute all the quantities needed for the gain
 % calculation
 
 gradStuff = interpOnTraj('unitGrad',testTraj,rayGd);
 neOnTraj = 10.^(interpOnTraj('valsNe',testTraj,rayGd));
 nc = rayBundleB1.nc(testIdx);
 neTonC = neOnTraj/nc;
 
 % compute vosc/c
 
 rayI15 = 1.e-15*rayBundleB1.I{testIdx};  % ray intensity in units
                                          % of 10^15 W/cm2
 omega0 = rayBundleB1.frequency(testIdx); % light frequency in rads/sec
 lam0 = cnst.twopi*cnst.c*1.e6/(rayBundleB1.frequency(testIdx)); % microns
 voscToc = (cnst.voscToC)*sqrt(rayI15)*lam0; % note shape is different!
 voscToc2 = voscToc.^2;                   % vosc^2/c^2
 
 % normalized scale length: (omega_0 L/c)
 %
 L = 1./gradStuff(:,3);                  % should check for zeros!
 Lnorm = (omega0*1.e-12/cnst.cumps)*L;   % dimensionless
 
 % density dependent part (alpha)
 %
 alphaFac1 = 1./(1.-sqrt(neTonC));  % see Rosenbluth formula in
                                    % Docs
 alphaFac2 = 1./sqrt(1-neTonC);   % ERROR!!!
 alphaFac = alphaFac1.*alphaFac2;
 
 otherFactors = alphaFac;          % need to add the things that
                                   % depend on k and theta/phi
 
 RosGain = (cnst.pi/4.0).*Lnorm.*(voscToc2').*otherFactors;
 
 
 figure(3)
 
 plot(rayBundleB1.trajs{testIdx}(:,1),rayBundleB1.I{testIdx})
 xlabel("time in ps")
 ylabel("intensity W/cm2")
 
 
 % inline function modified and moved to its own file (in ./Source)
 
 %function outval = dlnIdt(t,lnI,tSamp,gammaEM)
% outval = dlnIdt(t,lnI)
%   the temporal damping rate for EM waves
%
%    t     - time to evaluate the function
%    lnI   - unused, remove if ode45 is o.k. with it
%   tSamp  - the discrete times that gammaEM is known
%  gammaEM - the temporal intensity damping rate known at the
%            sampled times (tSamp)
%
%   returns the rhs to dlog(I)/dt = - gammaEM(t)
%
%  JFM: Sep 30, 2020
%

 %   outval = intper1(tSamp,gammaEM,t);

 %end
 