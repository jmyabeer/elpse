 function rayBundle = makeRayBundle(launchList,rayGd,segLength,itmax)

 % This function takes a launch list as an argument and returns a
 % ray bundle that can be integrated
 %    rayBundle is a cell array of ray structure
 %    launchList is a structure that contain
 %    initial conditions and other necessary information

 global cnst

 % [delta z, delta r] corresponding to 10 um (or segLength)
 % along the path
        
 if ~exist('segLength','var')
     segLength = 10;  % default (in um)
 end
 
 if ~exist('itmax','var')
     itmax = 5000;  % default
 end   
 
 
 % This function will initialize various types of beams based on
 % the "type"
 %
 bundleType = launchList.type;
 
 % We'll write this for 2-D then think how to generalize it to 3-D later
 %    add different cases as you need them...
 
 switch bundleType
   case 'laserBeam'
     nrays = launchList.nrays;
     omega = launchList.frequency*ones(1,nrays);     % rad/sec
     centroid = launchList.centroid;   % direction vector [kz,kr]
                                       % norm to 1
     focalPt = launchList.focalPt;
     spotShape = launchList.spot;      % just a diameter for now
                                        % (e.g. 500 um)
     % initial ray positions (col vecs) in beam wavefront coords
     rayWavPositions = zeros(2,nrays);  % each column is a ray
     rayWavPositions(2,:) = linspace(-spotShape.diameter/2, ...
                                     spotShape.diameter/2,nrays); % microns
     
     rayNorm = [-1,0];
     
     % can't do it this way (i.e. with dot); use atan instead
     %
     % Use switch to decide which quadrant you are in, then use
     % atan
     %
     angleB = atan2(centroid(2),centroid(1));
     angleA = atan2(0,-1);
     angle = angleB-angleA;

     % rotate to be parall to direction and translate to 5 mm away
     rotMat = [cos(angle),-sin(angle);sin(angle),cos(angle)];
     rotPos = rotMat*rayWavPositions;
     
     % translate origin to focus
     transPos = rotPos+repmat(focalPt',1,nrays);

     % translate 5 mm away (outside of the plasma)
     translate = -(launchList.translate)*centroid';   % microns
     rayPositions = transPos+repmat(translate,1,nrays);
     % transport into plasma can be done here?

     % set the frequency depending on the mode
     %
     switch launchList.mode
       case 'forward'
         rayBundle.frequency = omega;  % run time forwards   
       case 'backward'
         rayBundle.frequency = -omega; % run time backwards
       otherwise
         error('bad mode in launchList')
     end
     
     lambdaum = 1.e6*(cnst.twopi*cnst.c./omega);    % um
     rayBundle.nc = 1.1e21./lambdaum.^2;
     rayBundle.name = 'laserBeam';
     rayBundle.type = 'EM';
     rayBundle.mode = 'forward'; 

   case 'RamanAtDetector'
     nrays = launchList.nrays;
     omega = launchList.frequency*ones(1,nrays);     % rad/sec
     centroid = launchList.centroid;   % direction vector [kz,kr]
                                       % norm to 1
     focalPt = launchList.focalPt;
     spotShape = launchList.spot;      % just a diameter for now
                                        % (e.g. 500 um)
     % initial ray positions (col vecs) in beam wavefront coords
     rayWavPositions = zeros(2,nrays);  % each column is a ray
     rayWavPositions(2,:) = linspace(-spotShape.diameter/2, ...
                                     spotShape.diameter/2,nrays); % microns
     
     rayNorm = [1,0];
     angle = acos(dot(rayNorm,centroid)); % angle (rads) between beam and
                                          % simulation coords
     % rotate to be perp to direction and translate to 5 mm away
     rotMat = [cos(angle),-sin(angle);sin(angle),cos(angle)];
     rotPos = rotMat*rayWavPositions;
     
     % translate origin to focus
     transPos = rotPos+repmat(focalPt',1,nrays);

     % translate 5 mm away (outside of the plasma)
     translate = (launchList.translate)*centroid';   % microns
     rayPositions = transPos+repmat(translate,1,nrays);
     % transport into plasma can be done here?

     switch launchList.mode
       case 'forward'
         rayBundle.frequency = omega; % run time forwards   
       case 'backward'
         rayBundle.frequency = -omega; % run time backwards
       otherwise
         error('bad mode in launchList')
     end
     
     lambdaum = 1.e6*(cnst.twopi*cnst.c./omega);    % wavl in um
     rayBundle.nc = 1.1e21./lambdaum.^2;
     rayBundle.name = 'RamanAtDetector';
     rayBundle.type = 'EM';
     
   otherwise
     error('not a valid beam type')
 end


 rayBundle.nrays = nrays;
 rayBundle.direction = launchList.centroid;
 rayBundle.rayICs = rayPositions;
 rayBundle.trajs = cell(1,nrays);
 rayBundle.halt = zeros(1,nrays); % to terminate a trajectory

 % bring into the domain if required
 if strcmp(launchList.type,'laserBeam') | strcmp(launchList.type, ...
                                                 'RamanAtDetector')
     rayBundle = moveToDomain(rayBundle,rayGd,segLength,itmax);
 end

 end

 % helper subfunctions

 
function rayBundleOut = moveToDomain(rayBundle,rayGd,segLength,itmax)
% need to describe what this does
    
    % our quantum of path increment
    sgnFreq = sign(rayBundle.frequency(1));
    dispVec = sgnFreq*segLength*(rayBundle.direction)';
    
    inRange = getInRangeList(rayBundle,rayGd);
    its = 0;
    
    while ~all(inRange) & (its < itmax)
        % move all the out-of-domain rays 10 um
        dispVecRep = repmat(dispVec,1,sum(~inRange));
        colList = find(~inRange);
        rayBundle.rayICs(:,colList) = rayBundle.rayICs(:,colList) + ...
            dispVecRep;
        % get new list of out-of-domain rays
        inRange = getInRangeList(rayBundle,rayGd);
        its = its + 1;
    end

    % flag a potential problem (i.e., at least one ray never enters
    % the box)
    if its > itmax
        disp("exceeded max iterations in moveToDomain")
    end
    
    rayBundleOut = rayBundle;

end


function inRange = getInRangeList(rayBundle,rayGd)
% vectorized ...
    raysZ = rayBundle.rayICs(1,:);
    raysR = rayBundle.rayICs(2,:);
 
    inRrange = (raysR < rayGd.domain(4)) & (raysR > ...
                                            rayGd.domain(3));
    inZrange = (raysZ < rayGd.domain(2)) & (raysZ > ...
                                            rayGd.domain(1));
    inRange = inRrange & inZrange;
end

