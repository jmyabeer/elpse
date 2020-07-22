 % 
 % First cut at developing a ray integrator for EM and other waves:
 % (JFM June 17, 2020)
 %
 
 %
 % TO DO: figure out why the beam is smaller when you chose a
 % negative angle in the launch list. The differnce is large...
 %
 
 path(path,'/Plotting')
 path(path,'/Source')
 
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
 
 %
 % Get backscatter Raman on a chosen light trajectory
 %
 
 SRSangle = 336; % Raman scatter angle in degrees clockwise from the incident ray
 
 detector = initDetector;

 % Creates Raman for all rays in Beam 2 and checks if any Raman rays will
 % hit a detector
 detector.Source{1,1} = 'rayBundleB2';
 [srsBundle2,epwBundle2,detector] = Raman2Detector(rayBundleB2,rayGd,detector,SRSangle);

 % Creates Raman for all rays in Beam 1 and checks if any Raman rays will
 % hit a detector
 %detector.Source{2,1} = 'rayBundleB1';
 %[srsBundle1,epwBundle1,detector] = Raman2Detector(rayBundleB1,rayGd,detector,SRSangle);
 
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
 
