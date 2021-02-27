function plotOnRay(varToPlot,rayBundle,rayIdx,rayGd)
% function success = plotOnRay(varToPlot,rayBundle,rayIdx)
%       
%   This function plots a chosen 'hydro' variable along a single
%   ray belonging to a ray bundle
%
%      varToPlot - the 'hydro' variable: ''
%      rayBundle - a valid ray bundle containing the desired ray 
%         rayIdx - the index of the desired ray 
%          rayGd - the struct containing the hydro
%
%   Also: see interpOnTraj() 
%
% Last edited by: JFM 29/SEP/2020

global cnst

% Expand the following as we develop...
%
validVarList = ["valsNe","valsDLogNedz","valsDLogNedr","valsTe", ...
                "valsTi","gammaEM"];

% JFM: Try to be defensive - we need to revist because the variable
%  is not necessarily in rayGd (should check)
%  Some of the "derivative" ones are named differently - so we
%  should go back and fix that in interpOnTraj. While we're at
%  it we can include everyting - such as flow velocity and 
%  ionization
%
if ~contains(varToPlot,validVarList)
    error("Not a valid variable!")
end

if (rayIdx > rayBundle.nrays) | (rayIdx <= 0)
    error("Invalid ray index")
end


plotRay = rayBundle.trajs{rayIdx};   % the chosen trajectory
nCrit = rayBundle.nc(rayIdx);        % critical density for the ray

varOnRay  = interpOnTraj(varToPlot,plotRay,rayGd,nCrit);
timeOnRay = plotRay(:,1);            % time along ray in ps

% Plots

%hold on 

plot(timeOnRay,varOnRay)
xlabel("time in ps")
legend(varToPlot)

% need to finish this...

%hold off


end