 function varOnTraj = interpOnTraj(varToGet,y,rayGd,ncForRay)
 % Interpolates hydro variables onto a ray segment
 %   interpOnTraj(varToGet,y,rayGd,ncForRay)
 %             y  -   either cols constaing t z r kz kr
 %                    or cols containing just z r kz kr 
 %      varToGet  -  'valsNe', 'valsDLogNedz', 'valsDLogNedr'
 %                   'valsTe', 'valsTe', 'valsTi', 'gammaEM'
 %      ncForRay  -  critical density for the ray in cm^-3

    [~, ncols] = size(y);
    
    if ncols == 5 % this is the ode "traj" output - drop the time
        x = y(:,2:3);
    elseif ncols == 4 % no time
        x = y(:,1:2); % position vector is first two columns
    else
        error('wrong number of columns in y')
    end
    
    [ti,bc] = pointLocation(rayGd.DT,x);     % get Delauney
                                             % triangles
                                             
    % Thought - might add kLambadD to the below
    
    switch varToGet
      case 'valsNe'
        % to do - check to see if the variable exists
        triVals = rayGd.valsNe(rayGd.DT(ti,:));
      case 'valsDLogNedz'
        triVals = rayGd.valsDLogNedz(rayGd.DT(ti,:));
      case 'valsDLogNedr'
        triVals = rayGd.valsDLogNedr(rayGd.DT(ti,:));
      case 'valsTe'
        triVals = rayGd.valsTe(rayGd.DT(ti,:));
      case 'valsTi'
        triVals = rayGd.valsTi(rayGd.DT(ti,:));
      case 'gammaEM'
        % Electromagnetic intensity temporal damping rate
        %
        % electron temperature
        triValsTe = rayGd.valsTe(rayGd.DT(ti,:));
        TeOnTraj = dot(bc',triValsTe')';      % eV
        %
        % electron density
        triValsNe = rayGd.valsNe(rayGd.DT(ti,:));
        NeOnTraj = dot(bc',triValsNe')';
        NeOnTraj = 10.^NeOnTraj;              % cm^-3
        %
        % <Z>
        triValsZbar = rayGd.valsZbar(rayGd.DT(ti,:));
        ZbarOnTraj = dot(bc',triValsZbar')';
        %      
        % <Z^2>
        triValsZsqr = rayGd.valsZsqr(rayGd.DT(ti,:));
        ZsqrOnTraj = dot(bc',triValsZsqr')';
        %
        Z = ZsqrOnTraj./ZbarOnTraj;
        %
        % Let's assume that coulog = 8.0 (fix it later!)        
        coulog = 8.0;                  
        %
        nuei = (2.91e-6).*NeOnTraj.*coulog.*Z./(TeOnTraj.^(3/2)); % sec^-1
        %
        varOnTraj = (NeOnTraj./ncForRay).*nuei;  % sec^-1
      otherwise
        error("No matching hydro variable found in interpOnTraj");
    end
    
    if ~strcmp(varToGet,'gammaEM')
        varOnTraj = dot(bc',triVals')';
    end
    
 end
