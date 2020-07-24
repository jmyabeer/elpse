function haltFlags = checkPush(rayBundle,refracAllow)
% function haltFlags = checkPush(rayBundle,refractAllow)
%
%      rayBundle - a ray bundle (no assumption on type)
%   refractAllow - the amount in degrees in change between
%                  wavevectors for successive trajectory entries that triggers 
%                  the halt flag to be set (i.e. halt if the change
%                  is less than this.    
%
% Steven: checkPush determines whether each ray in the given ray bundle has been 
% pushed far enough so that the change in direction is within a given angle
% of allowance. If the ray's change in direction is within the allowance,
% a 1 will be put in the position in haltFlags corresponding to that ray.
% If the ray's change in direction is not within the allowance, a 0 will be
% put in the position in haltFlags corresponding to that ray.
% Currently will set halt flag if one of the changes in direction is more
% than the allowed change
%
%   JFM comment (24/July/2020): we should change this a bit. It
%   should reflect the amount of change per unit time (since the
%   successive trajectory entries could be very close together in
%   time).
%   Also - it should return a logical vector.
    

    % default angle of allowance
    if ~exist('refracAllow','var')
        refracAllow = 0.5; % degrees
    end
    
    for i =1:rayBundle.nrays
        % Calculates angle for last 5 positions in trajectory
        angle(:,i) = atand(rayBundle.trajs{1,i}(end-4:end,5)./ ...
            rayBundle.trajs{1,i}(end-4:end,4)); % atan(kr/kz)
    end
    
    % Calculates the change in angle between a position and the position
    % after
    %
    diff = angle(2:end,:) - angle(1:end-1,:);
    under = diff < refracAllow;
    
    % Sets halt flag if any change in direction is less than the allowed
    % change in direction
    %
    haltFlags = under(1,:)+under(2,:)+under(3,:)+under(4,:);
end
