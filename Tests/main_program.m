 % 
 % First cut at developing a ray integrator for EM and other waves:
 % (JFM June 17, 2020)
 %
 
 %
 % TO DO: figure out why the beam is smaller when you chose a
 % negative angle in the launch list. The differnce is large...
 %
 
 path(path,'./Plotting')
 path(path,'./Source')
 path(path,'./Steven_dev')
 
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

 % Beam 1
 %
 if ~exist('rayBundle','var')
     
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
     
 end

 % Beam 2
 %
 if ~exist('rayBundle','var')
     
     launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
     launchList.mode = 'forward';        % Could be backward also (neg omega?).
     launchList.nrays = 20;
     launchList.frequency = cnst.omega0; % 1/sec
     % center of spherical target
     launchList.focalPt = [-400,0];      % microns
     launchList.spot = struct('type','SG8','diameter',700); 
     angle = 180+(23.3); % (degres) is measured from "target norm"
     launchList.centroid = [cosd(angle),sind(angle)]; % unit vector in
                                                       % direction of
                                                       % beam propagation
                                                       % propagation
     launchList.translate = 5.0e3;   % distance in um from focus to
                                     % translate so that we are sure to
                                     % be far enough away to start
     % Create a ray bundle
     rayBundleB2 = makeRayBundle(launchList,rayGd);
     
 end

 %
 % Raman light bundle propagating back from a detector
 %

 if ~exist('ramanBundle1','var')
     ramanList1.type = 'RamanAtDetector';  
     ramanList1.mode = 'backward';   % integrate backwards (neg omega).
     ramanList1.nrays = 18;

     % derive frequency from detected wavelength in nm
     %
     detectWavl = (650.0)*1.e-9;    % meters
 
     ramanList1.frequency = cnst.twopi*cnst.c/detectWavl;  % 1/sec
     ramanList1.focalPt = [-400,0];                        % microns
     ramanList1.spot = struct('type','SG4','diameter',450);
     angle = 180-23.0;                                        % (degrees) 
     ramanList1.centroid = [cosd(angle),sind(angle)];      % unit vec
     ramanList1.translate = 5.0e3;   
     % Create a Raman ray bundle 
     %
     ramanBundle1 = makeRayBundle(ramanList1,rayGd);
 end

 if ~exist('ramanBundle2','var')
     ramanList2.type = 'RamanAtDetector';  
     ramanList2.mode = 'backward';   % integrate backwards (neg omega).
     ramanList2.nrays = 18;

     % derive frequency from detected wavelength in nm
     detectWavl = (700.0)*1.e-9;    % meters
 
     ramanList2.frequency = cnst.twopi*cnst.c/detectWavl;  % 1/sec
     ramanList2.focalPt = [-400,0];                        % microns
     ramanList2.spot = struct('type','SG4','diameter',450);
     angle = 180-23.0;                                        % (degrees) 
     ramanList2.centroid = [cosd(angle),sind(angle)];      % unit vec
     ramanList2.translate = 5.0e3;   
     % Create a Raman ray bundle 
     %
     ramanBundle2 = makeRayBundle(ramanList2,rayGd);
 end
 
 if ~exist('ramanBundle3','var')
     ramanList3.type = 'RamanAtDetector';  
     ramanList3.mode = 'backward';   % integrate backwards (neg omega).
     ramanList3.nrays = 18;

     % derive frequency from detected wavelength in nm
     detectWavl = (600.0)*1.e-9;    % meters
 
     ramanList3.frequency = cnst.twopi*cnst.c/detectWavl;  % 1/sec
     ramanList3.focalPt = [-400,0];                        % microns
     ramanList3.spot = struct('type','SG4','diameter',450);
     angle = 180+23.0;                                         % (degrees) 
     ramanList3.centroid = [cosd(angle),sind(angle)];      % unit vec
     ramanList3.translate = 5.0e3;   
     % Create a Raman ray bundle 
     %
     ramanBundle3 = makeRayBundle(ramanList3,rayGd);
 end
 
 %
 %  Integrate bundles
 %

 % incident laser light
 
 % First push
 %
 tPush = 3.7; % ps
 rayBundleB1 = pushBundle(rayBundleB1,rayGd,tPush,[100 100 100 100]);
 
 % refine
 %
 nits = 35;
 
 for i=1:nits
     rayBundleB1 = pushBundle(rayBundleB1,rayGd,0.1);
 end

  % incident laser light
 rayBundleB2 = pushBundle(rayBundleB2,rayGd,tPush,[100 100 100 100]);
 % refine
 %
 for i=1:nits
     rayBundleB2 = pushBundle(rayBundleB2,rayGd,0.1);
 end

 %
 % Halt any further integration of these rays
 %
 rayBundleB1.halt = setHaltAll(rayBundleB1);
 rayBundleB2.halt = setHaltAll(rayBundleB2);
 
 % Integrate the Raman light
 %
% $$$  ramanBundle1 = pushBundle(ramanBundle1,rayGd,4,[100 100 100 100]);
% $$$  % refine
% $$$  for i=1:20
% $$$      ramanBundle1 = pushBundle(ramanBundle1,rayGd,0.1);
% $$$  end

% $$$  ramanBundle2 = pushBundle(ramanBundle2,rayGd,4,[100 100 100 100]);
% $$$  % refine
% $$$  for i=1:20
% $$$      ramanBundle2 = pushBundle(ramanBundle2,rayGd,0.1);
% $$$  end
% $$$ 
% $$$  ramanBundle3 = pushBundle(ramanBundle3,rayGd,4,[100 100 100 100]);
% $$$  % refine
% $$$  for i=1:20
% $$$      ramanBundle3 = pushBundle(ramanBundle3,rayGd,0.1);
% $$$  end
 
 %
 % Get backscatter Raman on a chosen light trajectory
 %
 
 chosenRay = 13;   % pick a ray from the incident light bundle
 traj = rayBundleB2.trajs{chosenRay};
 freq = 1.e-12*rayBundleB2.frequency(chosenRay);  % ps^-1
 
 % TO DO: We should probably only choose a subset of the trajectory
 % points to integrate...
 
 % get backward SRS decay waves
 %
 [srsBundle,epwBundle] = getRamanWavevectors(traj,freq,rayGd);
  
 % Steven's version is called getRamanWavevector_sh()
 %  - it allows the backscatter angle to be chosen
 % 
 
 % TO DO: You should probably rename the above function to be called
 % "getRamanBundles()" or something like that
 
 % advance the Raman light
 %
 srsBundle2 = pushBundle(srsBundle,rayGd,2);
 
 % advance the Raman EPW
 %
 epwBundle2 = pushBundle(epwBundle,rayGd,10);
 
 
 %
 % Here you can try out Steven' "checkDetector()", and
 % "checkPush()" functions
 %
 
 %
 %  Make some plots
 %
 
 % Let's plot it on the density so that we can visualize the ode
 % solutions
 
 makePlotList(pltIncl,rayGd);

 % add the incident beam and SRS 
 %
 figure(1)
 hold on
 addBundlePlt(rayBundleB1,'b:');
 addBundlePlt(rayBundleB2,'b');
 %addBundlePlt(ramanBundle1,'r');
 %addBundlePlt(ramanBundle2,'g');
 %addBundlePlt(ramanBundle3,'m'); 

 % SRS light
 %
 rayRange = 1:12:srsBundle2.nrays;
 addBundlePlt(srsBundle2,'g',rayRange);
 
 % SRS EPW
 rayRange = 1:8:epwBundle2.nrays;
 addBundlePlt(epwBundle2,'r',rayRange);

 %
 %  Note to self: it looks like the EPW make a rainbow! 
 %    investigate
 %

 % Selected density contours
 %
 addContourPlt(rayBundleB1,rayGd,'nc');
 addContourPlt(rayBundleB1,rayGd,'nc4');
 addContourPlt(rayBundleB1,rayGd,'nc10');
 
 figure(2)
 hold on
 addBundlePlt(rayBundleB1,'b:');
 addBundlePlt(rayBundleB2,'b');
 
 % add a quiver plot of the flow velocity
 pltIncl = clearPlotList;
 pltIncl.quiverVel = true;
 makePlotList(pltIncl,rayGd);
 
