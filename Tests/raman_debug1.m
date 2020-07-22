 % 
 % Debugging the Raman routines
 %
 %   Edited: JFM 20/JUL/2020
 
 
 path(path,'./Plotting')
 path(path,'./Source')
 path(path,'./Steven_dev')
 
 global cnst
 
 % some initialization
 cnst = initCnst;             % will put more things in initCnst...
 pltIncl = initDefaultPlots;  % default plots to make
 pltIncl.temperature = true;  %   modify to add the elecron temp

  % get the hydro profile
 %
 dracoFile = "draco_EPsph_JFM.mat";

 % set the time slice here
 %
 tslice = 11;                

 % variables to import/define set here:
 %
 addVarFlag.ne = true;
 addVarFlag.dLogNedz = true;
 addVarFlag.dLogNedr = true;
 addVarFlag.Dmn = false;
 addVarFlag.te = true;
 addVarFlag.ti = false;
 addVarFlag.dLnTedz = true;
 addVarFlag.dLnTedr = true;
 addVarFlag.Vz = true;
 addVarFlag.Vr = true;
 
 % TO DO: add a function called "importAnalyticGrid()" that allows
 % you to define analytic hydro profiles, e.g. linearly varying
 % density, constant temperature and so on. Everything will still
 % be stored in the rayGd structure though.
 %
 
 if ~exist('rayGd','var')
     disp("loading hydro...")
     rayGd = importDracoGrid(dracoFile,tslice,addVarFlag);   
     disp("done loading hydro")
 else
     disp("using exisiting hydro")
 end
 
 % update if the time slice has changed too
 if rayGd.iTime ~= tslice
     disp("updating hydro...")
     rayGd = importDracoGrid(dracoFile,tslice,addVarFlag);
     disp("done updating hydro")
 end


 % Now just look at one problematic ray 
 
 SRSangle = 90;
 freq = 5.3663e+3;
 
 point = [5.2125 -65.0815 274.3752 -12.0888 9.2765];
 [srsBundle,epwBundle] = getRamanWavevectors_sh(point,freq,rayGd, ...
                                                      SRSangle);
 
 % set a breakpoint in getRamanWavevectors_sh and step through to
 % see where the complex numbers are coming from
 
 

 
 
 