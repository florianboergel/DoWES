%Function: vertical
%
%vertical takes a vector and gives back it in verical shape.
%--------------------------------
%USAGE:
%-------------
% [xv] = vertical(x)
%------------
%INPUT:
%-------------
% x     vector     in vertical or horizontal shape
%------------
%OUTPUT:
%-------------
% xv     vector     in vertical shape
%------------
%COMMENTS:
%-------------
% The code is replaced by the one of vertical_array
%------------
%Created by Juan JosE Trujillo on 
% On 2008-05-12-12:32
%----------------------------------
%
% Function: vertical_array.m
% It checks if the format of the array is vertical. If not it transpose it.
% --------------------------------
% Usage:
%-------------
% [ArrV]=vertical_array(Arr)
% [ArrV]=vertical_array(Arr,n_col)
% [ArrV,n_rows]=vertical_array(Arr,n_col)
%------------
% Input:
%-------------
% Arr    matrix  array
% ncol   float   number of expected columns (Default 1)
%------------
% Output:
%-------------
% ArrV   matrix  vertical array
% n_rows integer number of rows
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
% Davide Trabucchi om 2010-09-01
% It is now possible to use it as vertical.m
%------------
% Created: 
% Davide Trabucchi on 2010-08-18
% (c) Universitaet Oldenburg
% ---------------------------------

function [ArrV varargout]=vertical(Arr,varargin)

ncol=1;

if nargin==2
    ncol=varargin{1};
end

ArrDim=size(Arr);

test = ArrDim==ncol;

if test==[1 0]
    
    n_rows=ArrDim(2);
	ArrV=Arr';

elseif test
    n_rows=ArrDim(1);
    ArrV=Arr;
    warning('The proper format of the array could not be checked.')
    
elseif test==[0 1]
    
    ArrV=Arr;
    
else

    
	error(['The array has to be a Nx',num2str(ncol),' array.'])
	
    return
    
end

if nargout==2
    varargout{2}=n_row;
end

end




%Function: vertical
%
%vertical takes a vector and gives back it in verical shape.
%--------------------------------
%USAGE:
%-------------
% [xv] = vertical(x)
%------------
%INPUT:
%-------------
% x     vector     in vertical or horizontal shape
%------------
%OUTPUT:
%-------------
% xv     vector     in vertical shape
%------------
%COMMENTS:
%-------------
% 
%------------
%Created by Juan JosE Trujillo on 
% On 2008-05-12-12:32
%----------------------------------
% function [xv] = vertical(x)
% 
%   if isvector(x)
%     
%     if size(x,1) < size(x,2)
%       xv = x';
%     else 
%       xv = x;
%     end
%     
%   else
%      error('Input is not a vector!');
%   end
  