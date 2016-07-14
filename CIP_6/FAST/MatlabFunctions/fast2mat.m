% Function: fast2mat
% 
% This function reads in convert FAST .out/.elm output and loads it
% into matlab structure
%
% --------------------------------
% Usage:
%-------------
%[structFASTout] = fast2mat(fname)
%  
%------------
% Input:
%-------------
%
%fname           string      path and filename of FAST output file
%
%------------
%Output:
%-------------
%
%structFASTout   structure   
%
%------------
% Needs:
%-------------
%
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
% Marc Bromm on 2013-11-21-15:00
% (c) Universitaet Oldenburg
% ----------------------------------

function [structFASTout] = fast2mat(fname)

  %open file for input, include error handling
  fid = fopen(fname,'r');
  if fid < 0
    error(['Could not open ',fname,'.']);
  end
  
  cline = 1;
  tline = fgetl(fid);
  while ischar(tline);
   
    tlinens = regexprep(tline, '\s', '');

    %retrieve headers of data columns
    if(size(tlinens, 2) > 3)
      if(strcmp(tlinens(1:4), 'Time') == 1)    
        header_str = regexprep(tline, '\s', ' ');
        header = strsplit(header_str);      
        break;
      end
    end
    tline = fgetl(fid);    
    cline=cline+1;
  end
  
  %retrieve how many data columns there are
  coldata = size(header, 2);

  %read data
  tline = fgetl(fid);    
  data = fscanf(fid,'%f'); 
  datapoints = numel(data);
  
  if(mod(datapoints, coldata) ~= 0)
    error('unexpected data format')
  end

  rowdata = datapoints/coldata; 
  
  %arange data
  data = reshape(data, coldata, rowdata)';
  
  %build struct for output
  structFASTout = struct;
  structFASTout.header = header;
  
  for i = 1:coldata
     structFASTout.(header{i}) = data(:, i);
  end
 
end