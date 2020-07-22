 function varOnTraj = interpOnTraj(varToGet,y,rayGd)
 % Interpolates hydro variables onto a ray segment
 %   interpOnTraj(varToGet,y,rayGd)
 %             y  -   either cols constaing t z r kz kr
 %                    or cols containing just z r kz kr 
 %      varToGet  -  'valsNe', 'valsDLogNedz', 'valsDLogNedr'
 %                   'valsTe', 'valsTe', 'valsTi'

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
      otherwise
        error("No matching hydro variable found in interpOnTraj");
    end
    
    varOnTraj = dot(bc',triVals')';
 end
