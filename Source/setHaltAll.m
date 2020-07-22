function halt = setHaltAll(rayBundle,range)
% SETHALTALL - sets the halt flags to 1, or if range is specified
% only those with indices in the range are set e.g. range = [1 2 3
% 10];
%

    if ~exist('range','var')
        range = 1:rayBundle.nrays;
    end
    
    halt = rayBundle.halt;
    halt(range) = 1;
    
end
