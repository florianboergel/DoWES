function [WindData] = WindSlices(FileName, varargin)

% Function: WindSlices
% Function to visualize the magnitude of wind on different planes in a wind field. 
% --------------------------------
% Usage:
%-------------
% 1)
%    Windslices(FileName)
%    creates a plot with slices at all borders of the wind field, so-to-say
%    a wind field block
% 2)
%    Windslice(FilenName,'variable1', 'flag1/value1','variable2', 'flag2/value2')
%
%------------
% Input:
%-------------
% FileName      string   windfile in the workspace or path name of a
%                        windfile
% WindFlag      char     shown windspeed direction (longitudinal 'l',
%                        vertical 'v', horizontal 'h' and 'all')
% tslice        array    contains position of slicese at spec. time [s]
% yslice        array    contains position of slicese at spec. place [m]
% zslice        array    contains position of slicese at spec. place [m]
% tslice_num    int      number of periodic slices between the first and last slice
% tslice_i      array    contains position of slicese at spec. index
% hub           var      visualisation at hub center
%                        (1 = vertical plane, 2 = horizontal plane, 3 =
%                        both)
% tslice_bound  var      set or leave out boundary plane 
%                        (1 = front plane, 2 = last plane, 3 = both)
% trans         double   transparency of slices (1 = 100%, 0 = 0%)
% 
%------------
% Output:
%-------------
% "slice plot"
%
%------------
% local variables:
%-------------
% check           shows if a Wind File is loaded
% winddata        containing data of the Wind File
% v1,v            4-D wind velocity matrix 
% z,y,t           vectors, axes
% dz,dy,dt        distance between two points in the wind field [m]/[m]/[s]
% nz,ny,nt        number of points in the grid
% F               function to cal. the magnitude of wind at a spec. point
% offset          offset for hub height [m]
% ratio           visualisation ratio 
%
%------------
% Needs:
%-------------
% Wind File(.wnd)
% function: - readBLgrid
%           - VS2TS
%
%------------
% ToDo:
%-------------
%- write a polar coordinates optimized slice function
%------------
% Updated:
%-------------
% - works as well with INTs
%------------
% Created:
% Robert Menke 2012-12-11
%
% Updated:
% Bernd Kuhnle 2012-12-14
% Hauke Beck   2013-01-03
%
% (c) Universitaet Oldenburg
% ----------------------------------
% set Default
ratio      = 1.8;
hubheight  = 0;
ymeanstep  = 1;  % [m]  only for polar windfields
xmeanstep  = 1;  % [m]  only for polar windfields
tmeanstep  = 40; % [datasteps]  only for polar windfields

% ___read the varagins___
% count the varargins
optargin = size(varargin,2);
if optargin > 0
    for i = 1: optargin
        try
            if isequal(varargin{i}(end-2:end),'pas')
                pasfile  = varargin{i};
                varargin = varargin([1:i-1 i+1:length(varargin)]);   
            end
        end
    end               
end


% count the varargins
optargin = size(varargin,2);
% is the number of varagin odd or not
if rem(optargin,2) ~= 0
    error('You entered an odd number of inputs - Input has to entered pairwise')
end
% evaluate the varargins
if optargin > 0
    index_odd = 1 : 2 : optargin;   
    for i = 1: length(index_odd)
        eval([varargin{index_odd(i)}  ' = varargin{' num2str(index_odd(i)+1) '};']);
    end               
end

%___check whether a Wind File is in workspace or not___
try 
    WindData = evalin('caller', FileName);
    if isfield(WindData,'fc')
        filetype = 'wnd';
    elseif isfield(WindData,'nst')
        filetype = 'int';
    end
catch
    if exist(FileName,'file')
        if isequal(FileName(end-2:end),'wnd')
            filetype = 'wnd';   
            WindData = readBLgrid(FileName);
            % for Kaimal Windfields
            if WindData.fc == 7
                %  defrine the turbulence intensity for Kaimal
                if ~exist('tin','var')
                    tin = [10 8 5];
                end
                
                WindData = readBLgrid(FileName,'tin',tin);    
            end
        elseif isequal(FileName(end-2:end),'int')
            filetype = 'int';
            if exist(pasfile,'file')
                WindData = VS2TS(FileName, pasfile);
            else
                disp(['The file "' pasfile '" does not exist!' ])
            end
        end
    else
        disp(['The file "' FileName '" does not exist!' ])
        return
    end
end
    
switch filetype
    case 'wnd'    
        %___copy variables out of 'winddata'___
        v1 = WindData.velocity;
        nz = WindData.nz;
        ny = WindData.ny;
        nt = WindData.nt;
        dz = WindData.dz;
        dy = WindData.dy;
        dt = WindData.dx/WindData.MFFWS;
        zOffset = WindData.zOffset;

        %___change dimensions of the matrix___
        v = permute(v1,[3 1 4 2]); % [iy,it,iz,velocity]

        %___define the axes___
        z = 1:nz; 
        y = 1:ny;
        t = 1:nt;


        %___calculate the magnitude at a spec. point___
        if exist('WindFlag', 'var') == 0
            WindFlag = 'l'; % default flag
        end
        switch WindFlag
            case 'all'
                F = sqrt((v(y,t,z,1).^2)+(v(y,t,z,2).^2)+(v(y,t,z,3).^2));
            case 'l'
                F =v(y,t,z,1);
            case 'v'
                F =v(y,t,z,2);
            case 'h'
                F =v(y,t,z,3);
        end;
        
        %___offset, reference to the ground ___
        offset = zOffset-(max(z)*dz-dz)/2 ;

        %___axes ratio (y,z = 1)___
        %  ratio= (max(t)*dt-dt) / 360; % t

        %___set lateral zero to hub___
        y = ((-ny/2)+0.5):((ny/2)-0.5);
        
    case 'int'
 %___copy variables out of 'winddata'___
%         v1 = WindData.velocity;

        nt = WindData.n2t;

        dt = WindData.dtt;
        
        zOffset = 0;
        offset  = 0;
        
        v1 = WindData.data;
        
        % get coordinates of the grindpoints
        for i=1:WindData.nvt
            if i~=1
                theta(i) = floor((i-1)/WindData.nvt*WindData.nat)*(360/WindData.nat);
                index(i) = (i-(theta(i)/(360/WindData.nat)*WindData.nst)-1);
                rcalc(i) = WindData.rst(index(i));
            end
        end
        [int_y int_z] = pol2cartd(theta,rcalc);
        int_z = int_z + hubheight;
        
        int_yi = [nanmin(int_y):ymeanstep:nanmax(int_y)];
        int_zi = [nanmin(int_z):xmeanstep:nanmax(int_z)]';
        
        nz = length(int_zi);
        ny = length(int_yi);
        
        dz = xmeanstep;
        dy = ymeanstep;
        
        z = int_zi./dz; 
        y = int_yi./dy;
        t = 1:nt;

end

%___slices___

% yslice
if exist('yslice', 'var') == 0
    yslice = [];
else
    yslice = unique(yslice);
end

% zslice
if exist('zslice', 'var') == 0
    zslice = [];
else
    zslice = unique(zslice);
end

% tslice index (tslice_i)
if exist('tslice_i', 'var') == 0
    tslice_i = [];
else
    tslice_i = unique(tslice_i).*dt;
end

% tslice
if exist('tslice', 'var') == 0
    tslice = [];
else
    tslice = unique(tslice);
end

if exist('tslice_bound') == 1
    if tslice_bound == 1
        tslice = [tslice 0];
    elseif tslice_bound == 2
        tslice = [tslice max(t)*dt-dt];
    else
        tslice = [tslice 0 max(t)*dt-dt];
    end
    
end
        
if exist('tslice_num', 'var') == 1
    tslice_num =  (max(t)*dt-dt)/(tslice_num+1) : (max(t)*dt-dt)/...
        (tslice_num+1) : (max(t)*dt-dt)-(max(t)*dt-dt)/(tslice_num+1);
    tslice = [tslice tslice_num tslice_i];
else
    tslice = [tslice tslice_i];
end




%___hub___
if exist('hub') == 1
    if hub == 1
        yslice=[yslice 0];
    elseif hub == 2
        zslice=[zslice zOffset];
    else
        yslice=[yslice 0];
        zslice=[zslice zOffset];
    end
end

%__without variable input arguments___
if isempty(tslice) && isempty(yslice) && isempty(zslice)
    switch filetype
        case 'wnd'
            tslice = [tslice 0 max(t)*dt-dt];
            yslice = [yslice max(y)*dy min(y)*dy];
            zslice = [zslice offset offset+(max(z)-1).*dz];

        case 'int'
            tslice = [tslice 0 max(t)*dt-dt];
            yslice = [yslice max(y) min(y)];
            zslice = [zslice offset+(min(z)) offset+(max(z))];
            
    end

end
switch filetype
	case 'wnd'
        % check if the slices are well choosen
        tslice = tslice(tslice >= 0                       & tslice <= length(t)*dt);
        yslice = yslice(yslice >= min((y-1).*dy)          & yslice <= max((y-1).*dy));
        zslice = zslice(zslice >= min((z-1).*dz) + offset & zslice <= max((z-1).*dz) + offset);
        
    case 'int'
        % check if the slices are well choosen
        tslice = tslice(tslice >= 0                 & tslice <= max(t-1)*dt);
        yslice = yslice(yslice >= min((y))          & yslice <= max((y)));
        zslice = zslice(zslice >= min((z)) + offset & zslice <= max((z)) + offset);
end


switch filetype
    case 'int'
        if isempty(yslice) && isempty(zslice)
        	time4loop = unique(round((tslice./dt+1)));   
        else
            time4loop = unique(round([1 :tmeanstep :length(t) (tslice./dt+1) ]));
        end
        
        if ~isfield(WindData,'data_cart') 
                WindData.data_cart = ones(length(time4loop),1,length(int_yi),length(int_zi))*NaN;
                for i = 1:length(time4loop)
                    tic
                    WindData.data_cart(i,1,:,:) = griddata(int_y, int_z, v1(:,time4loop(i)), int_yi, int_zi,'cubic');
                    el(i) = toc;
                    if (WindData.n2t-i)*mean(el) > 30 && i == 20
                        disp(['Remaining time: ' sec2timestr((length(time4loop)-i)*mean(el))])
                    end
                end
        else
            if length(time4loop) ~= size(WindData.data_cart,1)
                WindData.data_cart = ones(length(time4loop),1,length(int_yi),length(int_zi))*NaN;
                for i = 1:length(time4loop)
                    tic
                    WindData.data_cart(i,1,:,:) = griddata(int_y, int_z, v1(:,time4loop(i)), int_yi, int_zi,'cubic');
                    el(i) = toc;
                    if (WindData.n2t-i)*mean(el) > 30 && i == 20
                        disp(['Remaining time: ' sec2timestr((length(time4loop)-i)*mean(el))])
                    end
                end
            end
        end
        t = time4loop;
        %___change dimensions of the matrix___
        F = permute(WindData.data_cart,[4 1 3 2]); % [iy,it,iz,velocity]

end

figure; % open plot in a new window
% create the slice-plot
switch filetype
    case 'wnd'
        SlicePlot = slice(t.*dt-dt,y.*dy,(z-1).*dz+offset,F,tslice,yslice,zslice);
    case 'int'
        SlicePlot = slice(time4loop.*dt-dt,y,z,F,tslice,yslice,zslice);
end


%___plot properties___
colormap jet; colorbar; 
xlabel('time [s]'); ylabel('width [m]'); zlabel('height [m]')
daspect([ratio 1 1]); % fits axes ratio
set(SlicePlot,'FaceColor','interp', 'EdgeColor','none')
switch filetype
    case 'wnd'
         switch WindFlag
            case 'all'
                componentstr = 'all';
            case 'l'
                componentstr= 'u';
            case 'v'
                componentstr= 'v';
            case 'h'
                componentstr= 'w';
        end;
        
        FileName_title = FileName;
        FileName_title(FileName=='_')= ' ';
        
        u_hh   = round(mean(F(ceil(WindData.ny/2),:,ceil(WindData.nz/2)))*100)/100;
        tin_hh = round(std(F(ceil(WindData.ny/2),:,ceil(WindData.nz/2)))/u_hh*100*100)/100;
        
        title({['W I N D F I E L D   V I S U A L I S AT I O N'],...
               [' '],...
               ['"' ' ' FileName_title(max(strfind(FileName_title, filesep))+1:end) ' ' '"'],...
               ['Component: ' componentstr '      U_h_u_b: ' num2str(u_hh) 'm/s      I_h_u_b: ' num2str(tin_hh) '%' ]});
           
    case 'int'
        componentstr = 'u';
        
        FileName_title = FileName;
        FileName_title(FileName=='_')= ' ';
        
        title({['W I N D F I E L D   V I S U A L I S AT I O N'],...
               [' '],...
               ['"' ' ' FileName_title(max(strfind(FileName_title, filesep))+1:end) ' ' '"'],...
               ['Component: ' componentstr '      U_h_u_b: ' num2str(round(WindData.QC.umean_hh*100)/100) 'm/s      I_h_u_b: ' num2str(round(WindData.QC.tin_hh*100*100)/100) '%' ],...
               ['U_W_F: '  num2str(round(WindData.QC.umean*100)/100) 'm/s      I_W_F: ' num2str(round(WindData.QC.tin*100*100)/100) '%' ] });
end


% limiting the plot
if isempty(yslice) && isempty(zslice)
    if length(unique(tslice)) == 1
    	xlim([min(tslice)-.1 max(tslice)+.1]) 
    else
        xlim([min(tslice) max(tslice)]) 
    end    
else
    xlim([0 max(t)*dt])
end

if isempty(tslice) && isempty(zslice)
    if length(unique(yslice)) == 1
    	ylim([min(yslice)-.1 max(yslice)+.1]) 
    else
        ylim([min(yslice) max(yslice)]) 
    end    
else
    ylim([min((y).*dy) max((y).*dy)])
end

if isempty(tslice) && isempty(yslice)
    if length(unique(zslice)) == 1
    	zlim([min(zslice)-.1 max(zslice)+.1]) 
    else
        zlim([min(zslice) max(zslice)]) 
    end
else
   zlim([min((z-1).*dz)+offset max((z-1).*dz)+offset]) 
end
    
% transparency
if exist('trans', 'var') == 1
    alpha(trans)
end



end

