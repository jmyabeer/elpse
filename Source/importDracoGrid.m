
 function rayGd = importDracoGrid(dracoFile,iTime,addVarFlag)
%
% Convert this Script into a function that returns a struct that
% holds the downsampled "blocks" i.e., zSblock, rSblock, the grid
% data, and the delauney triangulation for interpolation
%
% We will return:
%
%   rayGd - struct containing grid and data suitable for
%           interpolation
%
%
%    JFM April 15 (May 13), 2020
%
%
   
 ln10 = log(10); % don't want to rely on cnst being defined
 
 %dracoFile = "draco_N160421_001.mat";
 if ~exist('dracoFile','var')
     dracoFile = "draco_EPsph_JFM.mat";
 end
 
 % when this is a function addVarFlag should be passed as an argument
 if ~exist('addVarFlag','var')
     addVarFlag.ne = true;
     addVarFlag.dLogNedz = true;
     addVarFlag.dLogNedr = true;
     addVarFlag.Dmn = true;
     addVarFlag.te = true;
     addVarFlag.dLnTedz = true;
     addVarFlag.dLnTedr = true;
     addVarFlag.ti = true;
     addVarFlag.Vz = false;
     addVarFlag.Vr = false;
     addVarFlag.Zbar = true;
     addVarFlag.Zsqr = true;
 end
 
 
 % Load in the .mat file with hydro data
 T = load(dracoFile);
 
 times = T.Times;        % times are in ns 
   
 %
 % Data from all the time slices
 %

 % grid
 zAll = T.z;
 rAll = T.r;
 
 % electron density
 neAll = T.ne;
 % plasma temperatures
 teAll = T.Te;
 tiAll = T.Ti; 
 % velocity
 vzAll = T.Vz;
 vrAll = T.Vr;
 
 % average ionization <Z> and average square <Z^2>
 ZbarAll = T.Z;
 ZsqrAll = T.Zsq;
 
 
 % ------------ our time slice (check that it's legit)
 if iTime < 1
     iTime = 1
     disp("adjusting to first time slice")
 end
 if iTime > numel(times)
     iTime = numel(times)
     disp("adjusting to last slice")
 end

 rayGd.iTime = iTime;
 rayGd.time = times(iTime);
 
 % ------------
                  
 % use symmetry to complete grid for chosen time
 
 rSleft = -fliplr(rAll{iTime});
 rSAll = [rSleft,rAll{iTime}];        
 zSAll = [zAll{iTime},zAll{iTime}]; 
 
 % Subsampling:
 %   pick the stride to take over the original mesh (indices) for
 %   efficient subsampling. Might want to figure out how to
 %   choose/automate this sensibly

 [sizeInZ,sizeInR] = size(zAll{1});
 
 % choose sub-domain
 zstr = 1; 
 zstp = sizeInZ-1;
 rstr = 1;
 rstp = sizeInR;
 
 roff = 150/2;              % figure out how to fix automatically 
 zstp = round(0.8*zstp);    % OVERRIDE zstp index here
 
 % choose stride (could overwride this with optional func arg
 strdZ = 3;
 strdR = 3;
 
 % subsampled grid
 
 zSblock = zSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 rSblock = rSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);

 rayGd.zSblock = zSblock;
 rayGd.rSblock = rSblock;

 % [zmin zmax rmin rmax]
 rayGd.domain = [zSblock(1,1) zSblock(end,1) rSblock(1,1) ...
                   rSblock(1,end)];

 % subsample selected hydro variables for the chosen time

 if addVarFlag.ne
     % complete using symmetry for our time slice
     neSleft = fliplr(neAll{iTime});
     neSright = neAll{iTime}; 
     neSfull = [neSleft,neSright];
     % subsample
     neblock = log10(neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                         roff))-6; % log10 of density in 1/cm^3
     neblockUse = 1.e-6*neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                         roff); % 1/cm^3
     rayGd.neblock = neblock;
     % for interpolation
     rayGd.valsNe = reshape(neblock,[numel(zSblock),1]);
     
     % we might need the spatial derivaties too
     if addVarFlag.dLogNedz || addVarFlag.Dmn 
         [nZ nR] = size(neblock);
         dLogNedz = zeros(size(neblock));     % preallocate
                                              % d Log(density)/dz 
         for ridx = 1:nR
             % [dydx,~,~] = %
             % derivJFM(neblock(:,ridx),zSblock(:,ridx),'diff3');
             [dydx,~,~] = derivJFM(neblock(:,ridx),zSblock(:,ridx), ...
                               'pp');
             dLogNedz(:,ridx) = dydx;
         end
         % reshaped for interpolation
         rayGd.valsDLogNedz = reshape(dLogNedz,[numel(zSblock),1]);
     end
     
     if addVarFlag.dLogNedr || addVarFlag.Dmn
         dLogNedr = zeros(size(neblock));
         % d Log(density)/dr 
         for zidx = 1:nZ
             % [dydx,~,~] =
             % derivJFM(neblock(zidx,:),rSblock(zidx,:),'diff3');
             [dydx,~,~] = derivJFM(neblock(zidx,:),rSblock(zidx,:), ...
                               'pp');
             dLogNedr(zidx,:) = dydx;
         end
         % reshaped for interpolation
         rayGd.valsDLogNedr = reshape(dLogNedr,[numel(zSblock),1]);
     end
     
     % Once we have the derivative of ne we can compute the second
     % derivative and the cross term that are required for the
     % focusing tensor (JFM 27/MAY/2020)     
     if addVarFlag.Dmn
         d2Dzz = zeros(size(neblock));
         for ridx = 1:nR
             dnedz = neblockUse(:,ridx).*ln10.*dLogNedz(:,ridx);
             [dydx,~,~] = derivJFM(dnedz,zSblock(:,ridx),'pp');
             d2Dzz(:,ridx) = dydx;
         end
         % we probably need to do some smoothing/filtering on this
         % data
         rayGd.d2Dzz = d2Dzz;
         % reshaped for interpolation
         rayGd.valsD2Dzz = reshape(d2Dzz,[numel(zSblock),1]);
                  d2Dzz = zeros(size(neblock));
         
         d2Drr = zeros(size(neblock));        
         for zidx = 1:nZ
             dnedr = neblockUse(zidx,:).*ln10.*dLogNedr(zidx,:);
             [dydx,~,~] = derivJFM(dnedr,rSblock(zidx,:),'pp');
             d2Drr(zidx,:) = dydx;
         end
         % we probably need to do some smoothing/filtering on this
         % data
         rayGd.d2Drr = d2Drr;
         % reshaped for interpolation
         rayGd.valsD2Drr = reshape(d2Drr,[numel(zSblock),1]);

         d2Dzr = zeros(size(neblock));        
         for zidx = 1:nZ
             dnedz = neblockUse(zidx,:).*ln10.*dLogNedz(zidx,:);
             [dydx,~,~] = derivJFM(dnedz,rSblock(zidx,:),'pp');
             d2Dzr(zidx,:) = dydx;
         end
         % we probably need to do some smoothing/filtering on this
         % data
         rayGd.d2Dzr = d2Dzr;
         % reshaped for interpolation
         rayGd.valsD2Dzr = reshape(d2Dzr,[numel(zSblock),1]);      
     end

  end

 % electron temperature
 if addVarFlag.te
     teSleft = fliplr(teAll{iTime});
     teSright = teAll{iTime};
     teSfull = [teSleft,teSright];
     teblock = teSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.teblock = teblock;
     rayGd.valsTe = reshape(teblock,[numel(zSblock),1]);

     % we might need the spatial derivaties too
     if addVarFlag.dLnTedz 
         [nZ nR] = size(teblock);
         dLnTedz = zeros(size(teblock));     % preallocate
                                             % d Ln(Te)/dz 
         for ridx = 1:nR
             % [dydx,~,~] = %
             % derivJFM(teblock(:,ridx),zSblock(:,ridx),'diff3');
             [dydx,~,~] = derivJFM(log(teblock(:,ridx)),zSblock(:,ridx), ...
                               'pp');
             dLnTedz(:,ridx) = dydx;
         end
         % reshaped for interpolation
         rayGd.valsDLnTedz = reshape(dLnTedz,[numel(zSblock),1]);
     end

     if addVarFlag.dLnTedr 
         dLnTedr = zeros(size(teblock));
         % d Ln(Te)/dr 
         for zidx = 1:nZ
             % [dydx,~,~] =
             % derivJFM(teblock(zidx,:),rSblock(zidx,:),'diff3');
             [dydx,~,~] = derivJFM(log(teblock(zidx,:)),rSblock(zidx,:), ...
                               'pp');
             dLnTedr(zidx,:) = dydx;
         end
         % reshaped for interpolation
         rayGd.valsDLnTedr = reshape(dLnTedr,[numel(zSblock),1]);
      end
 end
 
 % ion temperature
 if addVarFlag.ti
     tiSleft = fliplr(tiAll{iTime});
     tiSright = tiAll{iTime};
     tiSfull = [tiSleft,tiSright];
     tiblock = tiSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.tiblock = tiblock;
     rayGd.valsTi = reshape(tiblock,[numel(zSblock),1]); 
 end

 % z-component of the flow velocity
 if addVarFlag.Vz 
     VzSleft = fliplr(vzAll{iTime});
     VzSright = vzAll{iTime};
     vzSfull = [VzSleft,VzSright];
     vzblock = vzSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.vzblock = vzblock;
     rayGd.valsVz = reshape(vzblock,[numel(zSblock),1]);      
 end

 % r-component of the flow velocity
 if addVarFlag.Vr
     VrSleft = fliplr(vrAll{iTime});
     VrSright = vrAll{iTime};
     vrSfull = [-VrSleft,VrSright]; % note the minus
     vrblock = vrSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.vrblock = vrblock;
     rayGd.valsVr = reshape(vrblock,[numel(zSblock),1]); 
 end
 
 % Zbar: average (over species) ionization state <Z>
 if addVarFlag.Zbar
     ZbarSleft = fliplr(ZbarAll{iTime});
     ZbarSright = ZbarAll{iTime};
     ZbarSfull = [ZbarSleft,ZbarSright]; 
     Zbarblock = ZbarSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.Zbarblock = Zbarblock;
     rayGd.valsZbar = reshape(Zbarblock,[numel(zSblock),1]); 
 end
 
 % Zsqr: average (over species) ionization state <Z^2>
 if addVarFlag.Zsqr
     ZsqrSleft = fliplr(ZsqrAll{iTime});
     ZsqrSright = ZsqrAll{iTime};
     ZsqrSfull = [ZsqrSleft,ZsqrSright]; 
     Zsqrblock = ZsqrSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
     rayGd.Zsqrblock = Zsqrblock;
     rayGd.valsZsqr = reshape(Zsqrblock,[numel(zSblock),1]); 
 end
 
   
 %
 %  Grid (points) for delaunay interpolation
 %
 
 PtsR = reshape(rSblock,[numel(rSblock),1]);
 PtsZ = reshape(zSblock,[numel(zSblock),1]);
 Pts = [PtsZ PtsR];

 % compute the triangulation
 rayGd.DT = delaunayTriangulation(Pts); 
 
 %  
 % Linear interpolation of data onto these points using a delaunay
 % triangulation query where Pq are the query points looks like 
 % 
 % [ti,bc] = pointLocation(DT,Pq);  
 % 
 % This finds the triangle that encloses each point - ti is the ID
 % and bc is barycentric coords
 %
 % The linear interpolation "VqLin" at the query point(s) Pq are
 % then given in the following way (using ValsNe as an example):
 %
 %   triVals = ValsNe(DT(ti,:));
 %   VqLin = dot(bc',triVals')';
 %

 end  % function

                                    
  
