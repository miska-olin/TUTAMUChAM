% Initializes the TUTAMUChAM
% TUTAMU Chamber Aerosol Model
% Miska Olin, August 2024

disp('TUTAMUChAM 1.0 starting...')
if ~exist('p','var')
    disp('Loading interpolation tables for the PL distribution ...')
    load interpolationTables
end

% Load the default setup
defaultSetupFile;

% Asking for own setup file
disp('Select the setup file')
filename = 'defaultSetupFile.m';
pathname = '';
[filename,pathname] = uigetfile('*.m','Select the setup file',strcat(pathname,filename));


% Runnig your own setup file
run(strcat(pathname,filename));

clear filename pathname

disp(p)




