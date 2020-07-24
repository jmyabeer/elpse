 % 
 % Debugging the Raman routines continued...
 %   Looks at a single ray from an incident bundle and
 %   makes Raman from it
 %       - use to test getRamanWavevectors_sh and the push
 %         of the generated SRS and Langmuir waves
 %   Edited: JFM 24/JUL/2020
 
 
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

  if ~exist('rayBundleB2','var')
     
     launchList.type = 'laserBeam';      % trigger for 'makeRayBundle'
     launchList.mode = 'forward';        % Could be backward also (neg omega?).
     launchList.nrays = 20;
     launchList.frequency = cnst.omega0; % 1/sec
     % center of spherical target
     launchList.focalPt = [-400,0];      % microns
     launchList.spot = struct('type','SG8','diameter',700); 
     angle = 180+(+23.3); % (degres) is measured from "target norm"
     launchList.centroid = [cosd(angle),sind(angle)]; % unit vector in
                                                       % direction of
                                                       % beam propagation
                                                       % propagation
     launchList.translate = 5.0e3;   % distance in um from focus to
                                     % translate so that we are sure to
                                     % be far enough away to start
     % Create a ray bundle
     rayBundleB2 = makeRayBundle(launchList,rayGd);

     % give it the right name
     rayBundleB1.name = 'Omega EP beam #2';     
 end

 
 %
 % Initialize the SABS detector
 %
 sabs = initDetector("SABS (Omega EP)");
 sabs.sourceBeams = [rayBundleB1];   % add all source beams before rays
                                     % are pushed
 sabs.sourceBeams = [sabs.sourceBeams rayBundleB2];
 
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
     rayBundleB1 = pushBundle(rayBundleB1,rayGd,0.1);
 end

 % Halt any further integration of these rays
 %
 rayBundleB1.halt = setHaltAll(rayBundleB1);

 
 %
 % Get backscatter Raman on a chosen light trajectory
 %
 
 % Define the current source parameters for our detector (time,
 % contributing bundle, SRS angle)
 %
 sabs.sourceParams{1,1} = rayGd.time;    % hydrodynamic time
 sabs.sourceParams{1,2} = 1;             % beam at index 1 in
                                         % sourceBeams array
 SRSangle = 220;  
 sabs.sourceParams{1,3} = SRSangle;      % Raman scatter angle in
                                         % degrees clockwise from the
                                         % incident ray
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
 addBundlePlt(rayBundleB1,'b:');

 
 % Loop over beam #1's rays
 %
 rSkip = 4;
 
 for chosenRay = 10:10 % 1:rSkip:rayBundleB1.nrays   

     traj = rayBundleB1.trajs{chosenRay};
     freq = 1.e-12*rayBundleB1.frequency(chosenRay);  % ps^-1
     
     % TO DO: We should probably only choose a subset of the trajectory
     % points to integrate...
     
     % get backward SRS decay waves
     
     [srsBundle,epwBundle] = getRamanWavevectors_sh(traj,freq,rayGd,SRSangle);
     
     % TO DO: You should probably rename the above function to be called
     % "getRamanBundles()" or something like that
     
     % advance the Raman light
     %
     srsBundle2 = pushBundle(srsBundle,rayGd,0.3);
     
     % advance the Raman EPW
     %
     epwBundle2 = pushBundle(epwBundle,rayGd,10);
     
     % Checks the change in trajectory of the SRS rays and further
     % integrates these rays if the change in direction is greater than the
     % allowed amount
     %
     while ~all(srsBundle2.halt)
         srsBundle2.halt = checkPush(srsBundle2);
         srsBundle2 = pushBundle(srsBundle2,rayGd,0.1);
     end
     
     % add SRS light rays to the plot
     %
     rayRange = 1:12:srsBundle2.nrays;
     addBundlePlt(srsBundle2,'g',rayRange);  

     % add SRS EPW to the plot
     rayRange = 1:12:epwBundle2.nrays;
     addBundlePlt(epwBundle2,'r',rayRange);

     % Score hits on the detector
     %  - think about how best to do this...
     sabs.frequencies{1,chosenRay} = checkDetector(srsBundle2,sabs);
     
 end
 

 % add some stuff to figure 2 for fun
 
 figure(2)
 hold on
 addBundlePlt(rayBundleB1,'b:');
 
 % add a quiver plot of the flow velocity
 pltIncl = clearPlotList;
 pltIncl.quiverVel = true;
 makePlotList(pltIncl,rayGd);
 
