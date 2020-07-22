%
%  Reads in Andrey's .mat file for EP spherical data
%    and takes derivatives of the hydro variables. 
%
%    JFM May 4, 2020
%
%  May 11 -> This is going to be the "init" function
%

 %dracoFile = "draco_N160421_001.mat";
 dracoFile = "draco_EPsph_JFM.mat";
 
 %saveFigs = true;
 saveFigs = false;
 
 % decide which plots to make
 plotFlags.qp = true;   % Quiver plot
 plotFlags.ne = true;
 plotFlags.te = false;
 plotFlags.ti = false;

 % figure out how many we are plotting
 ourFieldsAre = fieldnames(plotFlags);
 numToPlot = 0;
 for i = 1:length(ourFieldsAre)
     numToPlot = numToPlot+plotFlags.(ourFieldsAre{i});
 end
 currentFig = 1;
 
 % close all open figure windows
 close all
 
 % Load in the .mat file with hydro data
 
 if ~exist("T")
     T = load(dracoFile);
 end

 % need to have this
 %   put this stuff in an "initCnst" function (as globals)
 
 lambda0 = 0.351;             % microns
 nc = 1.1e21/lambda0^2;       % cm^-3
 loge10 = 2.302585092994046;  % natural log of 10
 
 %
 % Then we'll have "initHydro" with T global
 %
 
 times = T.Times;        % ns 
  
 
 % It looks like densities are in m^-3 - check with Andrey - YES
 
 % Most variables are stored in cell arrays, with a given time
 % corresponding to a single cell.
 %  - the contents of a cell are Matlab matrices
 
 
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
 

 % coordinate ranges to choose (need to figure out how to automate)
 
 [sizeInZ,sizeInR] = size(zAll{1});
 
 zstr = 1; 
 zstp = sizeInZ-1;
 rstr = 1;
 rstp = sizeInR;
 
 roff = 150/2;              % figure out how to fix automatically 
 zstp = round(0.8*zstp);    % OVERRIDE zstp index here
 
 
 % ------------ our time slice
 iTime = 8;      
 % ------------
                  
 
 % use symmetry to complete grid
 
 rSleft = -fliplr(rAll{iTime});
 rSAll = [rSleft,rAll{iTime}];        
 zSAll = [zAll{iTime},zAll{iTime}]; 

 % and the hydro variables for the chosen time

 if plotFlags.ne 
     neSleft = fliplr(neAll{iTime});
     neSright = neAll{iTime}; 
     neSfull = [neSleft,neSright];
 end
 
 if plotFlags.te
     teSleft = fliplr(teAll{iTime});
     teSright = teAll{iTime};
     teSfull = [teSleft,teSright];
 end
 
 if plotFlags.ti
     tiSleft = fliplr(tiAll{iTime});
     tiSright = tiAll{iTime};
     tiSfull = [tiSleft,tiSright];
 end
 
 
 % pick the stride to take over the original mesh (indices)

 %
 % "initBlock" will create the following block.zSblock etc.
 %
 
 strdZ = 3;
 strdR = 3;
 
 % subsampled grid
 
 zSblock = zSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 rSblock = rSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 
 % subsampled data
 
 if plotFlags.ne
     neblock = log10(neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                         roff))-6;
 end
 
 % also want this as the plasma frequency for interpolation
 
 if plotFlags.te
     teblock = teSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
 end

  % also want this as the plasma Debye length for interpolation
 
 if plotFlags.ti    
     tiblock = tiSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 end
 
  
 %
 % compute the z derivates
 %
 
 [nZ nR] = size(neblock);
 
 % preallocate
 dLogNedz = zeros(size(neblock));
 
 % d Log(density)/dz 
 for ridx = 1:nR
     % [dydx,~,~] =
     % derivJFM(neblock(:,ridx),zSblock(:,ridx),'diff3');
     [dydx,~,~] = derivJFM(neblock(:,ridx),zSblock(:,ridx),'pp');
     dLogNedz(:,ridx) = dydx;
 end

 %
 % compute the r derivates
 %
 
 % preallocate
 dLogNedr = zeros(size(neblock));
 
 % d Log(density)/dr 
 for zidx = 1:nZ
     % [dydx,~,~] =
     % derivJFM(neblock(zidx,:),rSblock(zidx,:),'diff3');
     [dydx,~,~] = derivJFM(neblock(zidx,:),rSblock(zidx,:),'pp');
     dLogNedr(zidx,:) = dydx;
 end
 
 % ----------------
 % we could interpolate onto a regular grid and make a quiver plot here
 % ----------------
 
 % Put in the grid points in the right form for interpolation
 PtsR = reshape(rSblock,[numel(rSblock),1]);
 PtsZ = reshape(zSblock,[numel(zSblock),1]);
 Pts = [PtsZ PtsR];

 % compute the triangulation of the grid points
 DT = delaunayTriangulation(Pts); 
 
 % Put the data in the right form for interpolation
 if exist("neblock")
     ValsNe = reshape(neblock,[numel(zSblock),1]); % Values on the
                                                   % points in 1-D array
 end
 % Derivatives
 if exist("dLogNedr")
     ValsDLogNedr = reshape(dLogNedr,[numel(zSblock),1]);
 end
 if exist("dLogNedz")
     ValsDLogNedz = reshape(dLogNedz,[numel(zSblock),1]);
 end
 
               
 % Setup query points [pointIdx,2] as an example for the
 % interpolation function
 %

 % linearly spaced grid
 %
 nskip = 50;     % microns
 rRange = [-rSblock(1,end):nskip:rSblock(1,end)];       % microns
 zRange = [zSblock(5,1):nskip:zSblock(end,1)];          % microns
 
 % Query points (points to interpolate onto)
 
 [ZQ,RQ] = meshgrid(zRange,rRange); % a test. We might want to do
                                    % logspace here and then use
                                    % this in computing the
                                    % derivatives
 

 % Query points as a linear array
 Pq = [reshape(ZQ,[numel(ZQ),1]) reshape(RQ,[numel(RQ),1])]; 
  
 % Linear interpolation of data onto these points using a Delaunay
 % triangulation query 
  
 [ti,bc] = pointLocation(DT,Pq);  % finds the triangle that
                                  % encloses each point - ti is the
                                  % ID and bc is barycentric coords
  
 % bring it back to a form where we can plot it
 if plotFlags.qp
     figure(currentFig)
     currentFig = currentFig+1;
     clf

     triVals1 = ValsDLogNedz(DT(ti,:));
     VqLin1 = dot(bc',triVals1')';

     triVals2 = ValsDLogNedr(DT(ti,:));
     VqLin2 = dot(bc',triVals2')';

     % plot the electron density
     pcolor(zSblock,rSblock,neblock)
     title(sprintf("Ne and velocity (arrows)",strdR,strdZ)) 
     xlabel("Z in microns")
     ylabel("R in microns") 
     shading interp  
     colorbar
 
     % plot velocity as arrows on top of density (vectors)
     hold on

     quiver(ZQ,RQ,-reshape(VqLin1,size(ZQ)),-reshape(VqLin2,size(ZQ)),'k') 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)),'Color','r', ...
          'FontSize', 14)
 
     hold off
    
     if saveFigs     
         saveas(gcf,"quiver-gradient","png")     
     end
 end

 
 
 % Then we will make a function "emRHS" that will do the
 % interpolation to create the RHS for the ray equations
 
 % Need to do temperature also
 

 if plotFlags.ne
     figure(currentFig)
     currentFig = currentFig+1;
     clf 

     % do a subplot here with density and its derivative
     hold on
     plot(zSblock(:,ridx),neblock(:,ridx),'b-')
     plot(zSblock(:,ridx),neblock(:,ridx),'bo')

     title(sprintf("Ne: lineout at R= %6.1f microns",rSblock(1,ridx))) 
     xlabel("Z in microns")
     ylabel("log10 n(x), dn/dx") 
 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)))
 
     if saveFigs     
         saveas(gcf,"density-deriv","png")     
     end
 end
 
 % temperatures
 
 if plotFlags.te      
     figure(currentFig)
     currentFig = currentFig+1;
     clf 

     plot(zSblock(ridx,:),teblock(ridx,:))

     title(sprintf("Te: lineout at R= %d",rSblock(1,ridx))) 
     xlabel("Z in microns") 
     ylabel("dTe/dx") 

     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)))
 
     if saveFigs    
         saveas(gcf,"temperature-elec-deriv","png")     
     end
 end
 

 if plotFlags.ti           
     figure(currentFig)
     currentFig = currentFig+1;
     clf 

     plot(zSblock(ridx,:),tiblock(ridx,:))

     title(sprintf("Ti: lineout at R= %d",rSblock(1,ridx))) 
     xlabel("Z in microns") 
     ylabel("dTi/dx") 

     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)))
 
     if saveFigs    
         saveas(gcf,"temperature-ion-deriv","png")     
     end 
 
 end

  
