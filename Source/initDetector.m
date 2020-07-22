function detector = initDetector(name,angPos,angAccept)
% func detetector = initDetector(name,angPos,angAccept) initializes 
%    the "detector" structure: 
%    
%      detector = struct(...         
%       'name', "default", ...       %  e.g. "FABS 32B" 
%       'angPos', [23,0], ...        %  [theta,phi] in degrees
%       'angAccept',3, ...           %  acceptance angle in degrees
%       'sourceBeams', [], ...       %  Array of rayBundle structs
%                                    %  corresponding to incident laser beams     
%       'sourceParams', [], ...      %  Each row describes the contribution to 
%                                    %  the corresponding cell of
%                                    %  frequecy field; labeled by time 
%                                    %  contributing bundle
%                                    %  (indexed into sourceBeams),
%                                    %  and SRS angle
% 
%       'histogram',[], ...          %  store histogram here
%                                    %  (reference by source)
%       'filter', []);               %  detection efficiency v. freq.
%       
%      detector.source{1,1} = 0.0;       % time in ns of hydro slice
%      detector.source{1,2} = "bundle";  % index of sourceBeams array
%      detector.source{1,3} = 180;       % SRS angle in degrees
    
    detector = struct(...         
     'name', "default", ...       %  e.g. "FABS 32B" 
     'angPos', [23,0], ...        %  [theta,phi] in degrees
     'angAccept',3, ...           %  acceptance angle in degrees
     'sourceBeams', [], ...       %  Array of rayBundle structs 
     'sourceParams', [], ...      %  Descriptor for each row of freqs      
     'frequencies',[], ...        %  store frequencies here
     'histogram',[], ...          %  store histogram here
     'filter', []);               %  detection efficiency v. freq.

    detector.sourceParams{1,1} = 0.0;       % time in ns of hydro slice
    detector.sourceParams{1,2} = 1;         % index of sourceBeams array
    detector.sourceParams{1,3} = 180;       % SRS angle in degrees     
    
    % Update name from arg
    if exist('name','var')
        detector.name = name;
    end        
        
    % Update angular position from arg
    if exist('angPos','var')
        detector.angPos = angPos;
    end
    
    % Update angular acceptance from arg
    if exist('angAccept','var')
        detector.angAccept = angPos;
    end
    
end