function DetectorFreq = checkDetector(srsBundle)

% checkDetector takes Raman trajectory and calculates the exit angle of
% each Raman ray to see if it will hit the detector (is within +/- 3
% degrees of the enterance angle of the incident beam [23 degrees] w.r.t 
% the z axis) and return a vector of frequencies of rays that will hit the
% detector to be used to create a histogram.

nrays = srsBundle.nrays;
j=1;

for i=1:nrays
    angle = atand(srsBundle.trajs{1,i}(end,5)/srsBundle.trajs{1,i}(end,4));
    if angle >= 20 && angle <=26
        DetectorFreq(j) = srsBundle.frequency(1,i);
        j = j+1;
    end
end

histogram(DetectorFreq)

end