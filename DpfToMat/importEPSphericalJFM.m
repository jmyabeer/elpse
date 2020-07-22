%
% Basic script to convert draco dpf files to .mat
%    From Andrey/Russ
%

if ~exist("getHydroOutputJFM")
    addpath(['/home/myatt/Jason-Things/ReseachTopics/NIF-' ...
             'experiments/MatlabDraco'],'-begin')
end    

workdir = '/home/myatt/Jason-Things/ReseachTopics/OMEGA-EP-experiments/';
shotdir = 'EPsph700um4nsRamp/';
matOutDir = shotdir;

hydroParams.shotNum=160421;
hydroParams.filename = strcat(workdir,shotdir,'draco.dpf');
hydroParams.readTimes = 'all';  %[6]; %'all'; %[6];   % Array of step indices or 'all'
hydroParams.runType='draco';
hydroParams.targetType='Planar';
hydroParams.CrossBeam=false;
hydroParams.TCC_offset = 293352;    % What is this?

%Filename to output to .mat file (optional)
matlabOutputFilename = strcat(workdir,matOutDir,'draco_EPsph_JFM.mat');


%% Import DPF file

T = getHydroOutputJFM(hydroParams.filename, hydroParams);
 
%% Save to .mat file 
%use T = load(<matlabOutputFilename>); to load

if exist('matlabOutputFilename', 'var')
    matlabOutputFilename = fixFilePath(matlabOutputFilename);
    save(matlabOutputFilename, '-struct', 'T');
end
