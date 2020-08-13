function DetectorFreq = checkDetector(rayBundle,Detector)

% checkDetector uses trajectory from rayBundle and calculates the exit angle of
% each ray to see if it will hit the detector (is within a given 
% acceptance of the enterance angle of the incident beam w.r.t the z axis) 
% and return a row vector of frequencies of rays that will hit the detector

% JFM comments (14/JUL/2020)
%   describe what DetectorFreq is (i.e. a row vector of frequencies
%   that correspond to the rays hitting the detector
   
    % Initializes hitIndxs row vector
    hitIndxs = [];     
    
    % Loop over rays in the rayBundle and collect the indices of
    % the rays in the rayBundle that hit the detector, storing them
    % in the array "hitIndxs":
    %
    %    - It might be useful to retrun hitIndxs
    %    
    for i=1:rayBundle.nrays
        % atan(kr/kz)
        %
        angle = atan2d(rayBundle.trajs{1,i}(end,5)/rayBundle.trajs{1,i}(end,4));
        if angle >= (Detector.angPos(1)-Detector.angAccept)...
                && angle <= (Detector.angPos(1)+Detector.angAccept)
            hitIndxs = [hitIndxs,i]; % Appends next hit ray index to hitIndxs
        end
    end
    
    % Checks if any rays hit the detector and returns the frequencies
    % corresponding rays that hit the detector
    %
    % JFM - find the freqencies of the rays indexed by "hitIndxs"
    %  and store them in the row vector "DetectorFreq"
    %
    if ~isempty(hitIndxs)
        DetectorFreq = rayBundle.frequency(hitIndxs);
    else
        DetectorFreq = [];
    end
end
