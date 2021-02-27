function outval = dlnIdt(t,lnI,tSamp,gammaEM)
% outval = dlnIdt(t,lnI,tSamp,gammaEM)
%   the temporal damping rate for EM waves
%
%    t     - time to evaluate the function
%    lnI   - unused, remove if ode45 is o.k. with it
%   tSamp  - the discrete times that gammaEM is known
%  gammaEM - the temporal intensity damping rate known at the
%            sampled times (tSamp)
%
%   returns the rhs to dlog(I)/dt = - gammaEM(t)
%
%  JFM: Sep 30, 2020
%

    gamma = -1.e-12*gammaEM;      % convert from s^-1 to ps^1 and
                                  % give the right sign
    
    if t < tSamp(1)
        outval = gamma(1);
    elseif t > tSamp(end)
        outval = gamma(end);
    else
        outval = interp1(tSamp,gamma,t);
    end
    
   % catch NaNs if any are encountered
   if isnan(outval)
       outval = 0;
   end
   
end