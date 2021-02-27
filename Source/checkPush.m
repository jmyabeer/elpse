function haltFlags = checkPush(rayBundle,refracAllow)
% function haltFlags = checkPush(rayBundle,refractAllow)
%
%      rayBundle - a ray bundle (no assumption on type)
%   refractAllow - the average amount in degrees of change per picosecond 
%                  between wavevectors for the last 11 entries below 
%                  which the halt flag will be set for the corresponding
%                  ray
%

    % default allowed rate of change in direction
    if ~exist('refracAllow','var')
        refracAllow = 0.5; % degrees/ps
    end
    
    for i =1:rayBundle.nrays
        % Calculates angle for last 11 positions in trajectory
        angle(:,i) = atan2d(rayBundle.trajs{1,i}(end-10:end,5), ...
            rayBundle.trajs{1,i}(end-10:end,4)); % atan(kr/kz)
        % Calculates the change in time (ps) between points for the last 5
        % trajectory points
        dt(:,i) = rayBundle.trajs{1,i}(end-9:end,1)- ...
                  rayBundle.trajs{1,i}(end-10:end-1,1);
    end
    
    % Calculates the change in angle between a position and the position
    % after
    dtheta = angle(2:end,:)-angle(1:end-1,:);
    
    % Calculates rate of change in direction between points
    dthetadt = dtheta./dt;
    
    % Sets halt flag to 1 if the average rate of change in direction is 
    % less than the allowed rate of change in direction
    haltFlags = mean(dthetadt) < refracAllow;
end
