% Script: Plot Turbulent Wind Fields
% Created to visualize turbulent wind fields created with TurbSim and used  
% in DOWES class (SS2014)
% by Luis, rev 0.0

% Matlab functions required:
% 'WindSlices','readBLgrid' & 'vertical' 

%% Main code

% add path with matlab functions
addpath('/media/fe/32492c60-8891-4a8f-8dc8-2a7888c7003a/Dokumente/DoWES/git_Dowes/CIP_5/CIP_Tutorial_5_-_Windfield_and_wake_simulation/TurbSim/MatlabFunctions');


% Inputs
timeinterval = [100,200,300,400,500,600]; % time sliced in timeseries [s]
yposition=0; %rotor center [m]
zposition=100; %rotor hub height [m]

%% Input turbulent file

% display dialog for file selection if no
[FileName, PathName] = uigetfile( 'TurbSim_test.wnd', 'Select full filed turbulent input file' );
% concatenate filename and path
InputFilePathAndName = strcat(PathName, FileName);

%% Plot 
% Call WindSlice with standard values to plot turbulent wind file
WindSlices(InputFilePathAndName,'tslice',timeinterval,'yslice',yposition,'zslice',zposition)
saveas(gcf, 'wind_wake_4_25ms.png')