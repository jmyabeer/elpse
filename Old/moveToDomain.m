
function rayBundleOut = moveToDomain(rayBundle,rayGd,segLength,itmax)
% need to describe what this does
    
    % [delta z, delta r] corresponding to 10 um (or segLength)
    % along the path
        
    if ~exist('segLength','var')
        segLength = 10;  % default (in um)
    end
    if ~exist('itmax','var')
        itmax = 5000;  % default
    end
    
    % our quantum of path increment
    sgnFreq = sign(rayBundle.frequency);
    dispVec = sgnFreq*segLength*(rayBundle.direction)';
    
    inRange = getInRangeList(rayBundle,rayGd);
    its = 0;
    
    while ~all(inRange) | (its > itmax)
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

% helper subfunction

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

