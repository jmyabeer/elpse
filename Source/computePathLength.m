 function pathLength = computePathLength(rayBundle,rayIdx)
 %   pathLength = computePathLength(rayBundle,rayIdx)
 %      rayBundle  -  a valid ray bundle struct
 %         rayIdx  -  (optional) returns the discrete path length
 %                    for this bundle, an index vector of bundles,
 %                    or all of them if absent
 %
 %  For example, use it like this:
 %
 %   >> rayBundle.path = computePathLength(rayBundle,[1:5])
 %
 %  JFM: Sep. 30, 2020.
 %

     nrays = rayBundle.nrays;
     
     if ~exist('rayIdx','var')
         rayIdx = [1:nrays];
     end
     
     if max(rayIdx) > nrays
         exit("invalid rayIdx")
     end
     
     pathLength = cell(1,nrays);           % preallocate 
     
     for ray = rayIdx
         traj = rayBundle.trajs{ray};      % get a ray trajectory
         [nr,~] = size(traj);
         currentPath = zeros(nr,1);        % preallocate column vector         
         dz2Anddr2 = diff(traj(:,2:3)).^2; % square of differences
         sumOfSqrs = sum(dz2Anddr2')';     % column vector
         currentPath(2:end,1) = cumsum(sqrt(sumOfSqrs));
         %
         pathLength{ray} = currentPath;    % package it up!
     end
     
 end
