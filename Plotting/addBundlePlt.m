function success = addBundlePlt(rayBundle,bcol,rayRange)
% ADDBUNDLETOPLT - to do later
%   

 if ~exist('rayRange','var')
     rayRange = 1:rayBundle.nrays;
 end
     
 for iplt = rayRange
     plot(rayBundle.trajs{iplt}(:,2),rayBundle.trajs{iplt}(:, ...
                                                           3),bcol)
 end
     
 success = 1;
 
end
