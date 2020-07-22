%
%  Reads in Andrey's .mat file obtained on using the Matlab script
%  "basicDracoImportJFM.m" and the appropriate Draco.dpf file.
%    This script make plots of density and electron temperature
%    as a function of time. Also creates the matlab
%    movie matrices Mcon (density) and Mcot (temperature). 
%
%    If you look at the struct "T" you will see all the variables
%    that are saved. There are many. We plot just ne and te here 
%    as an example.
%
%    JFM May 4, 2020
%
%

 dracoFile = "draco_EPsph_JFM.mat";
 
 %saveFigs = true;
 saveFigs = false;
 
 if ~exist("T")
     T = load(dracoFile);  % load the .mat file containing hydro
 end

 % laser parameters
 lambda0 = 0.351;        % microns
 nc = 1.1e21/lambda0^2;  % cm^-3
 qcrit = log10(nc/4);    % for plotting contour lines
 crit = log10(nc);       
 crit10 = log10(nc/10);  
 
 
 times = T.Times;        % Times corresponding to each slice in ns 
  
 
 % It looks like densities are in m^-3 - check with Andrey - YES
 
 % Most variables are stored in cell arrays, with a given time
 % corresponding to a single cell.
 %  - the contents of a cell are Matlab matrices
 
 
 % The coordinate grid is the same for all times (i.e., not Lagrangian)
 zAll = T.z;           
 rAll = T.r;
 neAll = T.ne;   % electron density time slices m^-3
 teAll = T.Te;   % electron temperature time slices (eV)
 
 % Although they are all the same, a coordinate grid is saved for
 % all time slices
 %
 [sizeInZ,sizeInR] = size(zAll{1}); % get the one for the first
                                    % time slice
 zstr = 1; 
 zstp = sizeInZ-1;
 rstr = 1;
 rstp = sizeInR;
 
 roff = 150/2;              % figure out how to fix automatically 

 zstp = round(0.8*zstp);    % OVERRIDE zstp index here

 rSleft = -fliplr(rAll{1});
 rSAll = [rSleft,rAll{1}];         
 zSAll = [zAll{1},zAll{1}];

 % pick a stride to take over the original mesh (indices)
 strdZ = 5;
 strdR = 5;
 
 % Take a subdomain, since the original one is a bit big
 %
 zSblock = zSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 rSblock = rSAll(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);

 
 %
 % Time dependent things follow
 %

 % density

 figure(1)
 clf 

 for iTime = 1:length(times)     
     % use symmetry to complete
     neSleft = fliplr(neAll{iTime});
     neSright = neAll{iTime};
     neSfull = [neSleft,neSright];  % append to make it look pretty...
     neblock = log10(neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2* ...
                             rstp-roff))-6;
 
     pcolor(zSblock,rSblock,neblock)
     hold on
     contour(zSblock,rSblock,neblock,[qcrit qcrit],'k:')
     contour(zSblock,rSblock,neblock,[crit crit],'k')  
     contour(zSblock,rSblock,neblock,[crit10 crit10],'k:')       
     hold off
     
     %axis([-450 800 -1000 1000])
     title(sprintf("Ne: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ))
     xlabel("Z in microns")
     ylabel("R in microns")
     shading interp 
     colorbar
     caxis([18 23])
 
     xt = 600; % microns
     yt = 1000;
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)))
     
     Mne(iTime) = getframe;
 end
 
 % to play, use: movie(Mne,nrepeats,fps) where fps is frames per second

 if saveFigs
     saveas(gcf,"density","png")
 end
 
 % temperature
 
 figure(2)
 clf 

 for iTime = 1:length(times)     
     % use symmetry to complete
     neSleft = fliplr(neAll{iTime});
     neSright = neAll{iTime};
     neSfull = [neSleft,neSright];  % append to make it look pretty...
     neblock = log10(neSfull(zstr:strdZ:zstp,rstr+roff:strdR:2* ...
                             rstp-roff))-6;     
     teSleft = fliplr(teAll{iTime});
     teSright = teAll{iTime};
     teSfull = [teSleft,teSright];
     teblock = teSfull(zstr:strdZ:zstp,rstr+roff:strdR:2*rstp-roff);
 
     pcolor(zSblock,rSblock,teblock)
     hold on
     contour(zSblock,rSblock,neblock,[qcrit qcrit],'k:')
     contour(zSblock,rSblock,neblock,[crit crit],'k')  
     contour(zSblock,rSblock,neblock,[crit10 crit10],'k:')       
     hold off
     
     %axis([-450 800 -1000 1000])
     title(sprintf("Te: original grid with stride %d in R and %d in Z", ...
               strdR,strdZ))
     xlabel("Z in microns")
     ylabel("R in microns")
     colormap hot
     caxis([200 2500])
     shading interp 
     colorbar
 
     xt = 600; % microns
     yt = 1000;
     text(xt,yt,sprintf("time: %3.2f ns",times(iTime)))

     Mte(iTime) = getframe;

 end
 
 % to play, use: movie(Mte,nreps,fps)  were fps is frames per
 % second
 