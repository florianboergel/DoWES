% Function: readBLgrid
% 
% "readBLgrid.m" reads .wnd files supported from GH-BLADED, AERODYN and FAST
% and return a structure with the windfield header and windfield data.
% "readBLgrid.m" supports Kaimal and v.Karman windfields
% --------------------------------
% Usage:
%-------------
% [ WindData ] = readBLgrid(filename, varargin)
%
% Example:
% for "von Karman" windfields
% [ WindData ] = readBLgrid('C:\Users\habe\WNDs\uref50_ti_4_s1.wnd')
% 
% for "Kaimal" windfields
% [ WindData ] = readBLgrid('C:\Users\habe\WNDs\uref50_ti_4_s1.wnd','ti',[5 4 2.5])
%                - scale the windfield to I_u = 5%
%                                         I_v = 4%
%                                         I_w = 2.5% turbulence intensity
%
% for "Mann" windfields
% [ WindData ] = readBLgrid('C:\Users\habe\WNDs\uref50_ti_4_s1.wnd','ti',[5])
%                - scale the windfield to I_u = 5%
%                                         I_v = is given by GH Bladed
%                                         I_w = is given by GH Bladed
%
%                                        
%------------
% Input:
%-------------
% filename      str             adress of wnd file
%------------
% Output:
%-------------
% WindData     struct           Header and windfielddata
%------------
% Needs:
%-------------
% - vertical.m
%------------
% ToDo:
%-------------
% 
%------------
% Updated:
%-------------
% - 2013-01-22: option to read in Kaimal windfields
% - 2013-02-13: option to read in Mann windfields
%------------
% Created:
% unkown - may be someone from NREL
% Changed:
% Bernd Kuhnle on 2013-02-13-11:46
% Hauke Beck on 2013-01-22-11:17
% (c) Universitaet Oldenburg
% ----------------------------------

function [WindData] = readBLgrid(FileName,varargin)
%% Read the VARARGINs
% Count the varargins
optargin = size(varargin,2);

% is the number of varagin odd or not 
if rem(optargin,2) ~= 0
    error('You entered an odd number of inputs - Input has to entered pairwise')
end
% Evaluate the varargins
if optargin > 0
        for i = 1 : optargin
            if isequal(varargin{i},'I') || isequal(varargin{i},'I0') || isequal(varargin{i},'i0') || isequal(varargin{i},'i') || isequal(varargin{i},'TI') || isequal(varargin{i},'ti')  || isequal(varargin{i},'Ti')...
                                        || isequal(varargin{i},'TIN') || isequal(varargin{i},'Tin') || isequal(varargin{i},'tin')
                tin = varargin{i+1};
            elseif isequal(varargin{i},'hh') || isequal(varargin{i},'HH') || isequal(varargin{i},'hubheight') || isequal(varargin{i},'Hubheight')
                hubheight = varargin{i+1};    
            end
        end
end
%% READ .WND BINARY FILE

len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end

% initialize variables
fileFmt  = 'int16';
ConvFact = 1.0; %results in meters and seconds

str      = {'HUB HEIGHT','CLOCKWISE','UBAR','TI(U','TI(V','TI(W'};  %MUST be in UPPER case
numVars  = length(str);   
SummVars = zeros(numVars, 1);

% open file
fid_wnd   = fopen( [ FileName '.wnd' ] );
if ( fid_wnd <= 0 )
   error( 'Wind file could not be opened.' );
   return;
end
 
nffc  = fread( fid_wnd, 1, 'int16' );                     % some kind of flag for BL & number of components
      
if nffc ~= -99  % AN OLD-STYLE AERODYN WIND FILE
    dz      = fread( fid_wnd, 1, 'int16' );               % delta z in mm
    dy      = fread( fid_wnd, 1, 'int16' );               % delta y in mm
    dx      = fread( fid_wnd, 1, 'int16' );               % delta x (actually t in this case) in mm
    nt      = fread( fid_wnd, 1, 'int16' );               % half number of time steps
    MFFWS   = fread( fid_wnd, 1, 'int16' );               % 10 times mean FF wind speed, should be equal to MWS
              fread( fid_wnd, 5, 'int16' );               % unnecessary lines
    nz      = fread( fid_wnd, 1, 'int16' );               % 1000 times number of points in vertical direction, max 32
    ny      = fread( fid_wnd, 1, 'int16' );               % 1000 times the number of points in horizontal direction, max 32
              fread( fid_wnd, 3*(-nffc-1), 'int16' );

   % convert the integers to real numbers 
    nffc     = -nffc;
    dz       = 0.001*ConvFact*dz;
    dy       = 0.001*ConvFact*dy;
    dx       = 0.001*ConvFact*dx;
    MFFWS    = 0.1*ConvFact*MFFWS;
    nz       = fix( mod(nz,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids
    ny       = fix( mod(ny,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids
        
else % THE NEWER-STYLE AERODYN WIND FILE
    fc       = fread( fid_wnd, 1, 'int16' );              % should be 4 to allow turbulence intensity to be stored in the header
    if fc == 4
        %% READING METHOD FOR v.KARMAN WINDFIELDS
        nffc         = fread( fid_wnd, 1, 'int32' );              % number of components (should be 3)
        lat          = fread( fid_wnd, 1, 'float32' );            % latitude (deg)
        z0           = fread( fid_wnd, 1, 'float32' );            % Roughness length (m)
        zOffset      = fread( fid_wnd, 1, 'float32' );            % Reference height (m) = Z(1) + GridHeight / 2.0
        TI_U         = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of u component (%)
        TI_V         = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of v component (%)
        TI_W         = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of w component (%)
        dz           = fread( fid_wnd, 1, 'float32' );            % delta z in m 
        dy           = fread( fid_wnd, 1, 'float32' );            % delta y in m
        dx           = fread( fid_wnd, 1, 'float32' );            % delta x in m           
        nt           = fread( fid_wnd, 1, 'int32' );              % half the number of time steps
        MFFWS        = fread( fid_wnd, 1, 'float32');             % mean full-field wind speed
        unused(1:3)  = fread( fid_wnd, 3, 'float32' );            % ??? unused variables (for BLADED)
        unused(4)    = fread( fid_wnd, 1, 'float32' );            % ??? unused variables (for BLADED)
        seed         = fread( fid_wnd, 1, 'int32' );              % Seed
        nz           = fread( fid_wnd, 1, 'int32' );              % number of points in vertical direction
        ny           = fread( fid_wnd, 1, 'int32' );              % number of points in horizontal direction
        for i = 1:3*(nffc-1)
            unused(end+1) = fread( fid_wnd, 1, 'float32' );       % ??? unused variables (for BLADED)
        end

        SummVars(3:6) = [MFFWS, TI_U, TI_V, TI_W];
       % set the amount of time steps
        nt       = max([nt*2,1]);
        % calculate the temporal resolution 
        dt       = dx/MFFWS;
        % the size of one time step
        nv       = nffc*ny*nz; 
        % calculate the scale factors
        Scale    = vertical(0.00001*SummVars(3)*SummVars(4:6))';
        % Windspeed offset in all directions
        Offset   = [SummVars(3) 0 0];
        % predefine some variables
        velocity = zeros(nt,nffc,ny,nz)*NaN;
        v_raw    = zeros(nffc, ny, nz, nt)*NaN;
        v        = nffc * ny * nz * nt;
        %% Reading the windfield
        % read the windfield
        [v cnt] = fread( fid_wnd, nv*nt, fileFmt );
         v_raw  = reshape(v,[nffc ny nz nt]);
         v_raw  = permute(v_raw,[4 1 2 3]);
        % scale the raw data
        velocity(:,:,:,:) = v_raw.*repmat(Scale,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]) + repmat(Offset,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]);
        %close the .wnd-file
        fclose(fid_wnd);


        WindData = struct('fc',fc,'nffc',nffc,'lat',lat,'z0',z0,'zOffset',zOffset,'TI_U',TI_U,...
        'TI_V',TI_V,'TI_W',TI_W,'dz',dz,'dy',dy,'dx',dx,'dt',dt,'nt',nt,'MFFWS',MFFWS,'nz',nz,'ny',ny,...
        'velocity',velocity,'Scale',Scale,'Offset',Offset,'seed',seed,'unused',unused);

    elseif fc == 7
        %% READING METHOD FOR KAIMAL WINDFIELDS
        % set the default Turbulence intensity [in %]
        if ~exist('tin','var')
            tin = [10 10 10].*[1 .8 .5];
        else
            if length(tin) ~= 3
                error('Please define the turbulence intensity for all components [Iu Iv Iw]')
            end
        end
        % set default hubheight [in m]
        if ~exist('hubheight','var')
            zOffset = 0;
        else
            zOffset = hubheight;
        end
        % set latitude
        if ~exist('latitude','var')
            lat = 53.143889; % [in °] - Oldenburg
        end
        % set roughness length
        if ~exist('z0','var')
            z0 = NaN; 
        end
        
        %% Reading the header
        unused(1)   = fread( fid_wnd, 1, 'int32' );                 % ??? unused variables (for BLADED)
        nffc        = fread( fid_wnd, 1, 'int32' );                 % the number of components (should be 3)
        dz          = fread( fid_wnd, 1, 'float32' );               % Resolution z in m
        dy          = fread( fid_wnd, 1, 'float32' );               % Resolution y in m
        dx          = fread( fid_wnd, 1, 'float32' );               % Resolution x in m
        nt          = fread( fid_wnd, 1, 'int32' );                 % Half the number of time steps
        MFFWS       = fread( fid_wnd, 1, 'float32' );               % Mean full-field wind speed in u
        unused(2)   = fread( fid_wnd, 1, 'float32' );               % ??? Mean full-field wind speed in v
        unused(3)   = fread( fid_wnd, 1, 'float32' );               % ??? Mean full-field wind speed in v
        xLu         = fread( fid_wnd, 1, 'float32' );               % Longitudinal turbulence Length scale
        unused(4)   = fread( fid_wnd, 1, 'float32' );               % ??? unused variables (for BLADED)
        seed        = fread( fid_wnd, 1, 'int32' );                 % Seed
        nz          = fread( fid_wnd, 1, 'int32' );                 % Number of Gridpoints in z
        ny          = fread( fid_wnd, 1, 'int32' );                 % Number of Gridpoints in y
        unused(5)   = fread( fid_wnd, 1, 'int32' );                 % ??? unused variables (for BLADED)
        unused(6)   = fread( fid_wnd, 1, 'float32' );               % ??? unused variables (for BLADED)
        xLv         = fread( fid_wnd, 1, 'float32' );               % Lateral turbulence Length scale
        unused(7)   = fread( fid_wnd, 1, 'float32' );               % ??? unused variables (for BLADED)
        unused(8)   = fread( fid_wnd, 1, 'float32' );               % ??? unused variables (for BLADED)
        xLw         = fread( fid_wnd, 1, 'float32' );               % Vertical turbulence Length scale
        COHDEC  	= fread( fid_wnd, 1, 'float32' );               % Coherence decay constat
        CohScale    = fread( fid_wnd, 1, 'float32' );               % Coherence Scale parameter [m]
        
        % set the amount of time steps
        nt       = max([nt*2,1]);
        % calculate the temporal resolution 
        dt       = dx/MFFWS;
        % the size of one time step
        nv       = nffc * ny * nz; 
        % Windspeed offset in all directions
        offset   = [MFFWS 0 0];
        % predefine some variables
        velocity = zeros(nt,nffc,ny,nz) * NaN;
        v_raw    = zeros(nffc, ny, nz, nt) * NaN;
        v        = nffc * ny * nz * nt;
        %% Reading the windfield
        % read the windfield
        [v cnt] = fread( fid_wnd, nv*nt, fileFmt );
         v_raw  = reshape(v,[nffc ny nz nt]);
         v_raw  = permute(v_raw,[4 1 2 3]);
        % calculate the real standart derivation
        std_real = mean(mean(std(v_raw),4),3);
        % calculate the scale factors
        Scale    = tin.*MFFWS./std_real./100;
        % scale the raw data
        velocity(:,:,:,:) = v_raw.*repmat(Scale,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]) + repmat(offset,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]);
        % close the .wnd-file
        fclose(fid_wnd);

        WindData = struct('fc',fc,'nffc',nffc,'MFFWS',MFFWS,'TI_U',tin(1),'TI_V',tin(2),'TI_W',tin(3),'zOffset',zOffset,'lat',lat,'z0',z0,...
        'dy',dy,'dz',dz,'dx',dx,'dt',dt,'nt',nt,'nz',nz,'ny',ny,...
        'velocity',velocity,'Scale',Scale,'Offset',offset,'seed',seed,'xLu',xLu,'xLv',xLv,'xLw',xLw,'CohDecay',COHDEC,'CohScale',CohScale,'unused',unused);

    elseif fc == 8
        %% READING METHOD FOR MANN WINDFIELDS
        if ~exist('tin','var')
            tin = [10].*[1];
        else
            if length(tin) ~= 1
                error('Please define the turbulence intensity for all components [Iu Iv Iw]')
            end
        end
        
        % Reading the header
        var(1)      = fread( fid_wnd, 1, 'int32' );           %1   
        nffc        = fread( fid_wnd, 1, 'int32' );       	    % nr. of components
        dz          = fread( fid_wnd, 1, 'float32' );           % delta z in m 
        dy          = fread( fid_wnd, 1, 'float32' );           % delta y in m
        dx          = fread( fid_wnd, 1, 'float32' );           % Alongwind spacing
        nt          = fread( fid_wnd, 1, 'int32' );             % half the number of time steps
        MFFWS       = fread( fid_wnd, 1, 'float32' );           % mean full-field wind speed
        var(end+1)  = fread( fid_wnd, 1, 'float32' );        %2
        var(end+1)  = fread( fid_wnd, 1, 'float32' );        %3
        var(end+1)  = fread( fid_wnd, 1, 'float32' );        %4   
        var(end+1)  = fread( fid_wnd, 1, 'float32' );        %5 
        seed        = fread( fid_wnd, 1, 'int32');              % mean full-field wind speed
        nz          =  fread( fid_wnd, 1, 'int32' );            % nr. of points in vertical direction
        ny          = fread( fid_wnd, 1, 'int32' );             % nr. of points in lateral direction
        var(end+1)  = fread( fid_wnd, 1, 'int32' );          %6
        var(end+1)  = fread( fid_wnd, 1, 'int32' );          %7 
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %8 
        var(end+1)  = fread( fid_wnd, 1, 'int32' );          %9 
        var(end+1)  = fread( fid_wnd, 1, 'int32' );          %10
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %11 
        gamma       = fread( fid_wnd, 1, 'float32' );           % shear 
        xLv         = fread( fid_wnd, 1, 'float32' );           % Longitudinal turbulence Length scale      
        xTI_V       =  fread( fid_wnd, 1, 'float32' );          % Scale for lateral turbulence          
        xTI_W       =  fread( fid_wnd, 1, 'float32' );          % Scale for vertical turbulence
        vert_FFT    =  fread( fid_wnd, 1, 'float32' );          % vertical size of blown up grid
        lat_FFT     =  fread( fid_wnd, 1, 'float32' );         	% lateral size of blown up grid
        nr_z_FFT    =  fread( fid_wnd, 1, 'int32' );            % nr of point in longitudinal direction for FFT
        nr_FFT      =  fread( fid_wnd, 1, 'int32' );            % nr. of points for FFT
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %13 
        var(end+1)  =  fread( fid_wnd, 1, 'float32' );       %14 
        var(end+1)  =  fread( fid_wnd, 1, 'float32' );       %15   
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %16   
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %17
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %18 
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %19 
        var(end+1)  =  fread( fid_wnd, 1, 'int32' );         %20  

       % set the amount of time steps
        nt       = max([nt*2,1]);
        % calculate the temporal resolution 
        dt       = dx/MFFWS;
        % the size of one time step
        nv       = nffc*ny*nz; 
        % calculate the scale factors
        Scale    = vertical(0.00001*MFFWS*[tin(1) tin(1).*xTI_V tin(1)*xTI_W])';
        % Windspeed offset in all directions
        Offset   = [MFFWS 0 0];
        % predefine some variables
        velocity = zeros(nt,nffc,ny,nz)*NaN;
        v_raw    = zeros(nffc, ny, nz, nt)*NaN;
        %% Reading the windfield
        % read the windfield
        [v cnt] = fread( fid_wnd, nv*nt+100, fileFmt );
        v_raw  = reshape(v,[nffc ny nz nt]);
        v_raw  = permute(v_raw,[4 1 2 3]);
        % scale the raw data
        velocity(:,:,:,:) = v_raw.*repmat(Scale,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]) + repmat(Offset,[size(v_raw,1) 1 size(v_raw,3) size(v_raw,4)]);
        %close the .wnd-file
        fclose(fid_wnd);
        
        WindData = struct('fc',fc,'nffc',nffc,'dz',dz,'dy',dy,'dx',dx,'dt',dt,'nt',nt,'MFFWS',MFFWS,'nz',nz,'ny',ny,...
        'velocity',velocity,'Scale',Scale,'Offset',Offset,'seed',seed,'gamma',gamma,'xLv',xLv,'xTI_V',xTI_V,'xTI_W',xTI_W,...
        'nr_FFT',nr_FFT,'vert_FFT',vert_FFT,'lat_FFT',lat_FFT,'nr_z_FFT',nr_z_FFT,'zOffset',0);
                
    end
end
return;