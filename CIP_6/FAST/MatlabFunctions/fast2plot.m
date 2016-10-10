% Function: fast2plot
% 
% This function reads FAST structure file "structFASTout" (output of 
% "fast2mat" function ) & plots each variable found in field "header" in a
% new created folder named "structFASTout-plot"
%
% --------------------------------
% Usage:
%-------------
% [messageFASTplot] = fast2plot(structFASTout)
%  
%------------
% Input:
%-------------
%
% structFASTout        string     structure including FAST outputs
%
%------------
%Output:
%-------------
%
% msgFASTplot          string     message indicating if plots are created
% plot files           jpg        plots in folder "structFASTout-plot"   
%
%------------
% Needs:
%-------------
% fast2mat
%------------
% ToDo:
%-------------
%
%------------
% Updated:
%-------------
%
%------------
% Created: 
% Luis Vera-Tudela on 2014-01-12-16:21
% (c) Universitaet Oldenburg
% ---------------------------------

function [messageFASTplot] = fast2plot(structFASTout)

%% Main Code
%clc, clear all

% number of variables to plot (first value is Time)
numVar = length(structFASTout.header)-1;

%numSteps = length(structFASTout.Time); % number of timesteps

%% Data handling

% folder name for plots   
folderName = 'Plots_TimeDomain_';
dateHour = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
figureFoldername = [folderName, dateHour];
    
% create folder if it doesn't exist  
if ~ (isdir(figureFoldername))
    mkdir(figureFoldername);
    cd(figureFoldername);
end

for fileNumber = 1:1:numVar
    
    % assign variable name (first one is Time)
    varName = ['structFASTout.', structFASTout.header{(1+fileNumber)}];
        
    % plot and save file per result
    
    % create figure
    figureHolder = figure('visible','off');
    % create axes
    axes1 = axes('Parent',figureHolder,'YGrid','on','XGrid','on');
    box(axes1,'on');
    hold(axes1,'all');
    % create plot vs Time
    plot(structFASTout.Time, eval(varName));
    % create labels
    xlabel('Time [s]');
    y_string = varName; % Place holder to accomodate units
    ylabel(y_string);
    % create file for varName
    plotFilename = [y_string, '.png'];
    saveas(figureHolder, plotFilename);
    close
end

% Return to folder with results
cd ../;

messageFASTplot = fprintf('Files created');
end