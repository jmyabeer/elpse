%
%  Reads in Andrey's .mat file for NIF shot N160421
%    make plots of density and electron temperature
%    as a function of time. 
%
%    JFM April 15 (May 4), 2020
%
%

 %dracoFile = "draco_N160421_001.mat";
 dracoFile = "draco_EPsph_JFM.mat";
 
 %saveFigs = true;
 saveFigs = false;
 makeInterpNePlot = false;
 makeQuiverPlot = true;
 
 % decide which plots to make
 plotFlags.ne = true;
 plotFlags.te = true;
 plotFlags.ti = false;
 plotFlags.Vz = false;
 plotFlags.Vr = false;

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
 
 lambda0 = 0.351;        % microns
 nc = 1.1e21/lambda0^2;  % cm^-3
 
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
 iTime = 16;      
 % ------------
                  
 
 % use symmetry to complete grid
 
 rSleft = -fliplr(rAll{iTime});
 rSAll = [rSleft,rAll{iTime}];        
 zSAll = [zAll{iTime},zAll{iTime}]; 

 % and the hydro variables for the chosen time

 if plotFlags.ne || makeQuiverPlot
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

 if plotFlags.Vz || makeQuiverPlot
     VzSleft = fliplr(vzAll{iTime});
     VzSright = vzAll{iTime};
     vzSfull = [VzSleft,VzSright];
 end

 if plotFlags.Vr || makeQuiverPlot
     VrSleft = fliplr(vrAll{iTime});
     VrSright = vrAll{iTime};
     vrSfull = [-VrSleft,VrSright]; % note the minus
 end
 
 
 % pick the stride to take over the original mesh (indices)

 strdZ = 5;
 strdR = 5;
 
 % subsampled grid
 
 zSblock = zSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 rSblock = rSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 
 % subsampled data
 
 if plotFlags.ne || makeQuiverPlot
     neblock = log10(neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                         roff))-6;
 end
 
 if plotFlags.te
     teblock = teSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp- ...
                       roff);
 end
 
 if plotFlags.ti    
     tiblock = tiSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 end
 
 if plotFlags.Vz || makeQuiverPlot
     vzblock = vzSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 end

 if plotFlags.Vr || makeQuiverPlot
     vrblock = vrSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 end
  
 %
 % Plot the original data
 %
 
 % density

 if plotFlags.ne
     figure(currentFig)
     currentFig = currentFig+1;
     clf 

     pcolor(zSblock,rSblock,neblock)
     %axis([-450 800 -1000 1000]) 
     title(sprintf("Ne: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ)) 
     xlabel("Z in microns")
     ylabel("R in microns") 
     shading interp  
     colorbar
 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)))
 
     if saveFigs     
         saveas(gcf,"density","png")     
     end
 end
 
 % temperatures
 
 if plotFlags.te      
     figure(currentFig)
     currentFig = currentFig+1;
     clf 

     pcolor(zSblock,rSblock,teblock) 
     %axis([-450 800 -1000 1000])
     title(sprintf("Te: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ)) 
     xlabel("Z in microns") 
     ylabel("R in microns") 
     colormap hot 
     shading interp  
     colorbar
 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)))
 
     if saveFigs    
         saveas(gcf,"temperature-elec","png")     
     end
 end
 

 if plotFlags.ti           
     figure(currentFig)
     currentFig = currentFig+1;
     clf 
 
     pcolor(zSblock,rSblock,tiblock) 
     %axis([-450 800 -1000 1000]) 
     title(sprintf("Ti: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ))
     xlabel("Z in microns")
     ylabel("R in microns") 
     colormap hot 
     shading interp  
     colorbar
 
     xt = 600; % microns 
     yt = 1000;
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)))

     if saveFigs    
         saveas(gcf,"temperature-ion","png")     
     end
 end

 % velocities
 
 if plotFlags.Vz           
     figure(currentFig)
     currentFig = currentFig+1;
     clf 
 
     pcolor(zSblock,rSblock,vzblock) 
     %axis([-450 800 -1000 1000]) 
     title(sprintf("Vz: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ))
     xlabel("Z in microns")
     ylabel("R in microns") 
     %     colormap hot 
     shading interp  
     colorbar
 
     xt = 600; % microns 
     yt = 1000;
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)))

     if saveFigs    
         saveas(gcf,"velocity-z","png")     
     end
 end

 
 if plotFlags.Vr           
     figure(currentFig)
     currentFig = currentFig+1;
     clf 
 
     pcolor(zSblock,rSblock,vrblock) 
     %axis([-450 800 -1000 1000]) 
     title(sprintf("Vr: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ))
     xlabel("Z in microns")
     ylabel("R in microns") 
     %     colormap hot 
     shading interp  
     colorbar
 
     xt = 600; % microns 
     yt = 1000;
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)))

     if saveFigs    
         saveas(gcf,"velocity-z","png")     
     end
 end
  
 %
 % should interpolate onto a regular grid to make a quiver plot
 %   quiver(zSblock,rSblock,vzblock,vrblock)
 %
 
 %
 %  Interpolation follows
 %
 
 % Put in the grid points in the right form for interpolation
 PtsR = reshape(rSblock,[numel(rSblock),1]);
 PtsZ = reshape(zSblock,[numel(zSblock),1]);
 Pts = [PtsZ PtsR];

 % compute the triangulation
 DT = delaunayTriangulation(Pts); 
 
 % Put the data in the right form for interpolation
 if exist("neblock")
     ValsNe = reshape(neblock,[numel(zSblock),1]); % Values on the
                                                   % points in 1-D array
 end

 if exist("vzblock")
     ValsVz = reshape(vzblock,[numel(zSblock),1]); % Values on the
                                                   % points in 1-D array
 end

 if exist("vrblock")
     ValsVr = reshape(vrblock,[numel(zSblock),1]); % Values on the
                                                   % points in 1-D array
 end
 
 % setup query points [pointIdx,2]
 %

 % linearly spaced grid
 %
 nskip = 80;     % microns
 rRange = [-rSblock(1,end):nskip:rSblock(1,end)];       % microns
 
% $$$  zLinMin = -400;         % microns
% $$$  zLinMax = 30 ;          % microns - must be positive!
% $$$  nZLin = 100;
% $$$  zLin = linspace(zLinMin,zLinMax,nZLin);   % linearly spaced portion
% $$$  deltaZ = zLin(end)-zLin(end-1);           % linear spacing
% $$$  
% $$$  nZLog = 40;
% $$$  zLogMax = 800;  % end of logarithmic domain in z direction (microns)
% $$$  zLogRange = logspace(log10(zLinMax+deltaZ),log10(zLogMax),nZLog);
 
% $$$  zRange = [zLin zLogRange];   % concatenate linear and log grids

 zRange = [zSblock(5,1):nskip:zSblock(end,1)];   % microns
 
 % Query points (points to interpolate onto)
 
 [ZQ,RQ] = meshgrid(zRange,rRange); % a test. We might want to do
                                    % logspace here and then use
                                    % this in computing the
                                    % derivatives
 

 Pq = [reshape(ZQ,[numel(ZQ),1]) reshape(RQ,[numel(RQ),1])]; 
  
 % Linear interpolation of data onto these points using a delaunay
 % triangulation query 
  
 [ti,bc] = pointLocation(DT,Pq);  % finds the triangle that
                                  % encloses each point - ti is the
                                  % ID and bc is barycentric coords
  
 % bring it back to a form where we can plot it
 if makeInterpNePlot
     figure(numToPlot+1)
     clf

     triVals = ValsNe(DT(ti,:));
     VqLin = dot(bc',triVals')';

     % plot the density on the new grid to see how it looks
     pcolor(ZQ,RQ,reshape(VqLin,size(ZQ))) 
     %     shading interp
     title("density on new grid")
     xlabel("Z in microns")
     ylabel("R in microns")
     colorbar
 end

 % This is probably the more useful one
 if makeQuiverPlot
     figure(numToPlot+1)
     clf

     triVals1 = ValsVz(DT(ti,:));
     VqLin1 = dot(bc',triVals1')';

     triVals2 = ValsVr(DT(ti,:));
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

     quiver(ZQ,RQ,reshape(VqLin1,size(ZQ)),reshape(VqLin2,size(ZQ)),'k') 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ps",times(iTime)),'Color','r', ...
          'FontSize', 14)
 
     hold off
    
     if saveFigs     
         saveas(gcf,"quiver-vels","png")     
     end
 end
  
