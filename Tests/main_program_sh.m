 % 
 % First cut at developing a ray integrator for EM and other waves:
 % (JFM June 17, 2020)
 %   Edited: JFM 14/JUL/2020
 
 
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

 % incident laser light - beam 1
 
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

 % incident laser light - beam 2
 %
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
 
 %SRSangle = 180; % Raman scatter angle in degrees clockwise from the incident ray
 SRSangle = 150;  % Raman scatter angle in degrees clockwise from
                  % the incident ray
 
 % Debugging notes: run into problems with interpolation at +140 degrees
 %
 
 detector.AngPos    = 23; % Degrees couterclockwise w.r.t positive z axis 
 detector.AngAccept = 3;  % Degrees
 
 detector.Source{1} = 'rayBundleB2';
 
 % Loop over beam #2's rays
 %
 rSkip = 4;
 
 for chosenRay = 1:rSkip:rayBundleB2.nrays   
    traj = rayBundleB2.trajs{chosenRay};
    freq = 1.e-12*rayBundleB2.frequency(chosenRay);  % ps^-1
 
    % TO DO: We should probably only choose a subset of the trajectory
    % points to integrate...
 
    % get backward SRS decay waves
 
    [srsBundle,epwBundle] = getRamanWavevectors_sh(traj,freq,rayGd,SRSangle);
  
    % TO DO: You should probably rename the above function to be called
    %   "getRamanBundles()" or something like that
    %
    %    Also, need to remove the Raman from the outbound light
    %    trajectory
    %
 
    % advance the Raman light
    
    srsBundle2 = pushBundle(srsBundle,rayGd,2.0);
 
    % advance the Raman EPW
    
    epwBundle2 = pushBundle(epwBundle,rayGd,10.0);
    
    % Checks the change in trajectory of the SRS rays and further
    % integrates these rays if the change in direction is greater than the
    % allowed amount
    %
    while ~all(srsBundle2.halt)
        % set halt flag if outbound rays are good
        %
        srsBundle2.halt = checkPush(srsBundle2);
        % ... and push the ones that need it
        %
        srsBundle2 = pushBundle(srsBundle2,rayGd,0.5); % think
                                                       % about push time
    end
    
    % checkDetector returns a row vector of frequencies of rays in
    % the srsBundle2 that hit the detector (or an empty vector)
    %
    detector.Freq{1,chosenRay} = checkDetector(srsBundle2,detector);
    
 end
 
 % Now do the same thing for the laser beam #1
 %
 detector.Source{2,1} = 'rayBundleB1';
 
 rSkip = 4;
 
 for chosenRay = 1:rSkip:rayBundleB1.nrays   
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
    srsBundle2 = pushBundle(srsBundle,rayGd,2);
 
    % advance the Raman EPW
    %
    epwBundle2 = pushBundle(epwBundle,rayGd,10);
    
    % Checks the change in trajectory of the SRS rays and further
    % integrates these rays if the change in direction is greater than the
    % allowed amount
    %
    while ~all(srsBundle2.halt)
      srsBundle2.halt = checkPush(srsBundle2);
      srsBundle2 = pushBundle(srsBundle2,rayGd,0.5);
    end
    
    detector.Freq{2,chosenRay} = checkDetector(srsBundle2,detector);
    
 end
 
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
 
