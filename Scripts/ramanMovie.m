 % 
 % Makes a movie showing SRS backscatter one ray at a time
 % (JFM June 18, 2020)
 %
 
 path(path,'./Plotting')
 path(path,'./Source')
 
 global cnst
 
 % some initialization
 cnst = initCnst;             % will put more things in initCnst...
 pltIncl = initDefaultPlots;  % default plots to make
 pltIncl.temperature = false;  %   modify to add the elecron temp
 
 
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
 % initialize and create a launch list and ray bundles 
 %

 % Beam 1
 %
 if ~exist('rayBundleB1','var')
     
     launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
     launchList.mode = 'forward';        % Could be backward also (neg omega?).
     launchList.nrays = 20;
     launchList.frequency = cnst.omega0; % 1/sec
     % center of spherical target
     launchList.focalPt = [-400,0];      % microns
     launchList.spot = struct('type','SG8','diameter',750); 
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
 if ~exist('rayBundleB2','var')
     
     launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
     launchList.mode = 'forward';        % Could be backward also (neg omega?).
     launchList.nrays = 20;
     launchList.frequency = cnst.omega0; % 1/sec
     % center of spherical target
     launchList.focalPt = [-400,0];      % microns
     launchList.spot = struct('type','SG8','diameter',750); 
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
 %  Integrate bundles
 %

 % incident laser light
 
 tPush = 3.7; % ps
 rayBundleB1 = pushBundle(rayBundleB1,rayGd,tPush,[100 100 100 100]);
 rayBundleB2 = pushBundle(rayBundleB2,rayGd,tPush,[100 100 100 ...
                     100]);
 % refine
 %
 nits = 35;
 
 for i=1:nits
    rayBundleB1 = pushBundle(rayBundleB1,rayGd,0.1);
    rayBundleB2 = pushBundle(rayBundleB2,rayGd,0.1);     
 end

 % Halt any further integration of these rays
 %
 rayBundleB1.halt = setHaltAll(rayBundleB1);
 rayBundleB2.halt = setHaltAll(rayBundleB2);
 
 %
 % Get backscatter Raman on a chosen light trajectory
 %
 
 for chosenRay = 1:rayBundleB2.nrays
     % chosenRay = 13;   % pick a ray from the incident light
     % bundle
     %
     traj = rayBundleB2.trajs{chosenRay};
     freq = 1.e-12*rayBundleB2.frequency(chosenRay);  % ps^-1
 
     % get backward SRS decay waves
     %
     [srsBundle,epwBundle] = getRamanWavevectors(traj,freq,rayGd);
   
     % advance the Raman light and EPW
     %
     srsBundle2 = pushBundle(srsBundle,rayGd,2);
     epwBundle2 = pushBundle(epwBundle,rayGd,10);
     
     %  Make some plots
     %      Let's plot it on the density 
 
     makePlotList(pltIncl,rayGd);

     % add the incident beam and SRS 
     %
     hold on

     addBundlePlt(rayBundleB1,'b:');
     addBundlePlt(rayBundleB2,'b');

     % SRS light
     %
     rayRange = 1:12:srsBundle2.nrays;
     addBundlePlt(srsBundle2,'g',rayRange);
 
     % SRS EPW
     rayRange = 1:8:epwBundle2.nrays;
     addBundlePlt(epwBundle2,'r',rayRange);

     addContourPlt(rayBundleB1,rayGd,'nc');
     addContourPlt(rayBundleB1,rayGd,'nc4');
     addContourPlt(rayBundleB1,rayGd,'nc10');
     
     Mraman(chosenRay) = getframe;
 end
 
 % play movie with: movie(Mraman,nreps,fps)
 
