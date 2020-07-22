function DetectorFreq = checkDetector_jfm(srsBundle)

% checkDetector takes Raman trajectory and calculates the exit angle of
% each Raman ray to see if it will hit the detector (is within +/- 3
% degrees of the enterance angle of the incident beam [23 degrees] w.r.t 
% the z axis) and return a vector of frequencies of rays that will hit the
% detector to be used to create a histogram.

nrays = srsBundle.nrays;

%
%  Ideas: I think we need a "preprocessor" which grabs the last several k-vectors for each
%  trajectory and tests to see how much they are changing. It will
%  then return the indices of the rays that need to be pushed more with
%  pushBundle. We can set the halt flag on those that don't need to
%  be moved. Only when they're all good do we call this
%  function. What do you think?
%
%  It might be useful to have a "detector" score card. This could
%  be a struct which records the index of rays that have hit
%  the detector and the name of the contributing bundle (we should add a name
%  field to the rayBundle structs). Remember that we'll end up
%  having a bundle for a range of SRS backscatter angles, so we'll
%  end up with a list of bundle names that contribute. The
%  DectectorFreq array that you define below could then be a field of this
%  detector - perhaps belonging to a cell array (one for each 
%  contributing bundle. Thinking ahead - we will need to do this for every
%  time that we have hydro data for, so it should also identify the
%  time (for when we make a streaked detector image). I suggest we
%  have a different detector score card for each hydro profile
%  (time) otherwise things are going to get (more) horribly complicated.

%
%  If you pass the detector score card as an argument, then it
%  could also contain its angle and angular acceptance. That way
%  you won't need to hard code these things (and it will be easier
%  to change them or use them for other detectors like those on the
%  NIF).
%
%  Also, we might want to separate out the plotting from the
%  computation. By this, I mean would could write some separate
%  plotting functions that take the detector score card as input
%  and make the appropriate plots.
%

%  Have a think about the above and let's discuss on Monday (June
%  22).

j=1;

for i=1:nrays
    angle = atand(srsBundle.trajs{1,i}(end,5)/srsBundle.trajs{1,i}(end,4));
    if angle >= 20 && angle <=26
        DetectorFreq(j) = srsBundle.frequency(1,i);
        j = j+1;
    end
end

% how about vectorizing this and/or just returning the indices of
% the rays that hit? This could be a two-liner:

%angle = atand(srsBundle.trajs{1,:}(end,5)/srsBundle.trajs{1,:}(end,4));
%hitIndxs = find( (angle >= 20) && (angle <= 26));

% the frequencies are then given as needed by
%DetectorFreq = srsBundle.frequency(hitIdxs);
% although you will only want to do this if the index list is non-empty

histogram(DetectorFreq)

end