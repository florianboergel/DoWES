function [file_path_and_name] = select_file_from_folder(PathName)
% Function:SELECT_FILE_FROM_FOLDER
% It allows interactively to select a file from a given (default) folder
% --------------------------------
% Keywords: file, folder, selection
%-------------
% Usage:
% [file_path_and_name] = select_file_from_folder(PathName)
%
% Examples:
% 
% myfile = select_file_from_folder(absolute_path)
%------------
% Input:
%-------------
% PathName - absolute path to select file from
%------------
% Output:
%-------------
% file_path_and_name - file path and name as variable
%------------
% Dependencies:
%-------------
% none
%------------
% ToDo:
%-------------
% 
%------------
% Updated:
%-------------
% Luis Vera-Tudela on 2013-11-05
%------------
% Created: 
% Luis Vera-Tudela on 2013-06-23
% (c) Universitaet Oldenburg
% ----------------------------------

%% Main code

% display dialog for file selection
[FileName, PathName] = uigetfile(PathName,'Select input file');

% concatenate filename and path
file_path_and_name = strcat(PathName, FileName);

% Check if file exist
CheckFileExist = exist(file_path_and_name,'file'); % (Should be equal to 2)

if ~CheckFileExist
    error(['File ',file_path_and_name,' was not found!']);
else
    disp(['File ',file_path_and_name,' was found']);
end