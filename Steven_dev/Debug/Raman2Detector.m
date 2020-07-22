% Not sure what this function should be called
function [srsBundle2,epwBundle2,detector] = Raman2Detector(rayBundle,rayGd,...
                                                        detector,SRSangle)
% Loops through all rays in "rayBundle" and creates Raman scattering
% events. Raman light is then integrated far enough so that the
% diffraction due to change in density is small enough to get an accurate
% exit angle for all Raman rays. These exit angles are then used to 
% determine if a ray will hit the detector and stores the frequencies of
% rays that hit in a row vector.

    persistent counter
    
    if isempty(counter)
        counter = 1;
    else
        counter = counter+1;
    end
    
%     for chosenRay = 1:rayBundle.nrays   
    chosenRay = 13;
        traj = rayBundle.trajs{chosenRay};
        freq = 1.e-12*rayBundle.frequency(chosenRay);  % ps^-1
 
        % get SRS decay waves
 
        [srsBundle,epwBundle] = getRamanBundles(traj,freq,rayGd,SRSangle);
  
        % TO DO: Remove the Raman from the outbound light trajectory
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
    [detector.Freq{counter,chosenRay}, detector.hitIndxs{counter,chosenRay}]...
        = checkDetector(srsBundle2,detector);
    
%     end % for loop through all rays
end