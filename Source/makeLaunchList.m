 
 function launchList = makeLaunchList(launchType)
 % This function takes a launch type  as an argument and returns a
 % launch list that can be used to create ray initial conditions
 %    rayBundle is a cell array of ray structure
 %    launchList is a cell array of launch structures that contain
 %    initial conditions and other necessary information

 % We'll have functions that create a launch list for two-plasmon
 % decay, the various forms of SRS, and LDI (i.e., a launch list for
 % Langmuir waves).

 % launch lists can also be created if we encounter mode conversion
 % points for both EM and LWs.

 
 launchList = cell(100,1);
 
 switch launchType
   case 'laser beam'
     % a drive beam (do this one first)
     % need direction, spot size, and poining
     % also frequency
   case 'langmuir'
     % launch some plasma waves
   case 'ramanLight'
     % we have some SRS to trace
   case 'twoPlasmon'
     % this might be a useful special case
   otherwise
     error("launchType not valid in makeLaunchList")
 end

 
 end   % function


 
 