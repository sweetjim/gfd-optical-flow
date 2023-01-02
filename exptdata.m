classdef exptdata<handle
    %EXPTDATA is a class that stores information and data about the
    %loaded experiment and to conduct velocity field calculations using
    %Optical Flow algorithms.

    properties  % Global variables 
        folder  % Folder information about loaded NetCDF file
        ncinfo  % NetCDF file information
        image   % Last image output
        dimensions = struct('ds',1/5e-3,'dt',1) % Space (ds,m) and time (dt,s) interval dimensions
        times   % Time array
        velocity % Optical-Flow velocity outputs
    end

    methods % CONSTRUCTOR
        function obj = exptdata(folderHandle)
            % EXPTDATA Constructor
            % Upon construction, the following public properties will be
            % initialized:
            %   folder, ncinfo, times
            %_________________________________________________________________
            %_________________________________________________________________
            %   Input argument:
            %   folderHandle    -   [string OR DIR struct]
            %                       Path designation to .nc file or DIR
            %                       struct of same property.
            %_________________________________________________________________
            %_________________________________________________________________
            %   Output argument:
            %   obj             -   [dynamic class]
            %                       This class object.
            %_________________________________________________________________
            %_________________________________________________________________
            %   Examples:
            %   ~~~ 1 ~~~ 
            %   file = uigetfile;       %   will return a pathname string
            %   >> file = 'example_folder\folder_with_data\example_file.nc'
            %
            %   data = exptdata(file);  %   will assign 'data' to the
            %                               workspace
            %   ~~~ 2 ~~~
            %   files    = uigetdir;
            %   files    = dir(files,'*.nc');
            %   >> files =  [nx1 struct array with fields]
            %   >>          name
            %   >>          folder
            %   >>          date
            %   >>          bytes
            %   >>          isdir
            %   >>          datenum
            %
            %   data    = exptdata(files(i)) %  where i is an index between
            %                                   [0,n].
            %_________________________________________________________________
            %_________________________________________________________________
            %%
            if isa(folderHandle,'char')
                folderHandle = dir(folderHandle);   % conversion to DIR struct
            end
            obj.folder  = folderHandle;
            obj.ncinfo  = ncinfo(fullfile(obj.folder.folder,obj.folder.name));
            obj.times   = seconds([1:obj.ncinfo.Dimensions(1).Length]-1);
        end
    end
    methods % VELOCITY
        function plotVelocity(obj,imageindex,varargin)
            %% PLOTVELOCITY is a tool for viewing velocity outputs
            % NOTE: After parsing this function, a variable termed
            % "velocity" will appear in the base-workspace.
            %_________________________________________________________________
            %_________________________________________________________________
            % Input arguments:
            %   REQUIRED:
            %   imageindex  -  [int] index of the time-series to be loaded
            %   
            %   OPTIONAL:
            %   Argument parsing: (argval)
            %
            %   'mask'      -   [bool]
            %                   Masking state;
            %                   true OR 1 to apply mask [DEFAULT]
            %                   false OR 0 otherwise
            %
            %   OPTIONAL:
            %   Argument parsing: (argname,argval)
            %   ARGNAME     |   ARGVAL
            %   ------------|-------------
            %   'ax'        -   [axes handle] 
            %                   DEFAULT IS GCA (GET CURRENT AXES)
            %
            %   'input'     -   [struct] 
            %                   Velocity structure outputted from
            %                   "getVelocity" function.
            %                   DEFAULT IS LAST OUTPUT (IN PROPERTIES)
            %
            %   'view'      -   [string] 
            %                   Viewing mode selection:
            %                   'Image'; raw input [DEFAULT]
            %                   'u'; horz. velocity field
            %                   'v'; vert. velocity field
            %                   'U'; velocity mag. field
            %                   'Vorticity'; vorticity field
            %
            %   'overlay'   -   [string]
            %                   Overlay mode selection:
            %                   'None'; No overlay [DEFAULT]
            %                   'Vectors'; Quiver overlay
            %                   'Streamlines'; Streamline overlay
            %
            % 
            %   'clim'      -   [1x2 vec, double]
            %                   C-axis limits. 
            %                   DEFAULTS ARE VARIABLE 
            %                   (FOUND UNDER 'Viewing modes' and 'Labels'
            %                   sections):
            %                       [25 30] - 'Image' mode 
            %                       [0 inf] - 'U' mode
            %                       [-1 1].*maxval - otherwise
            %_________________________________________________________________
            %_________________________________________________________________
            %   Example:
            %   % Preamble
            %   data = exptdata(filename);
            %   imageindex = 100;
            %   fig = figure;
            %   ax  = axes(fig);
            %   data.getVelocity(imageindex,'method','slow');
            %   
            %   data.plotVelocity(imageindex,'ax',ax,'u','Vectors')
            %   % this will show the unmasked u-velocity output with a
            %   % quiver field overlay.
            %_________________________________________________________________
            %_________________________________________________________________
            %% Input parsing
            ax      = [];               % Axes to plot to
            in      = obj.velocity;     % Velocity structure
            view    = 'Image';          % Viewing mode
            overlay = 'None';           % Overlay mode
            mask    = true;             % Masking state
            parseInput(varargin)
            
            if isempty(ax)
                ax = gca;
            end

            if isempty(in)
                error('No velocity field has been calculated.')
            end
            %% Parsing
            im      = in.im1;
            [X,Z,t,ds] = obj.getDimensions;

            ux  = in.u;
            uy  = in.v;
            [zs,xs] = size(ux); % If resizing has occured, size(X)~=size(ux,1) etc.
            Z = linspace(-max(Z)/2,max(Z)/2,zs); % Symmetric vectors
            X = linspace(-max(X)/2,max(X)/2,xs);
            
            if mask
                ux = obj.applyMask(ux);
                uy = obj.applyMask(uy);
            end

            U   = sqrt(ux.^2+uy.^2);
            vor = vorticity(ux,uy);
            %% Viewing modes
            adata = mask;
            cmap = 'balance';
            quivercolor = 'r';
            switch view
                case 'Image'
                    img = im;
                    adata = mask;
                    cmap = 'parula';
                    quivercolor = 'k';
                    clim = [25 30];
                case 'u'
                    img = ux;
                case 'v'
                    img = uy;
                case 'U'
                    img = U;
                    clim = [0 inf];
                    cmap = 'gray';
                case 'Vorticity'
                    img = vor;
            end
            
            imagesc(ax,X,Z,img,'AlphaData',adata)
            
            axis(ax,'tight','xy');
            %% Overlay modes
            XLIM = xlim(ax);
            YLIM = ylim(ax);
            hold(ax,'on')
            switch overlay
                case 'Vectors'
                    h = vis_flow(ux,uy,...
                            'points',40,...
                            'offset',1, ...
                            'mag',3,...
                            'col','m',...
                            'ax',ax,...
                            'X',X,...
                            'Z',Z);
                        set(h, 'Color', quivercolor);
                case 'Streamlines'
                    h=streamslice(ax,X,Z,ux,uy,2);
                        set(h, 'Color', 'r');
            end            
            %% Labels
            xlim(ax,XLIM)
            ylim(ax,YLIM)
            xlabel(ax,'x (m)');
            ylabel(ax,'z (m)');
            
            
            t0          = obj.times(1);
            addlabels('ax',ax,'x','x (m)')
            title(ax,sprintf(...
                    't = %s (hh:mm:ss)',...
                    datestr(diff([t0 t(imageindex)]),'HH:MM:SS')));
            hold(ax,'off')
            
            if strcmp(cmap,'balance')
                clim = max(img,[],'all').*[-1 1];
            end
            caxis(ax,clim)
            addColorbar('ax',ax,'cmap',cmap,'pivot',0)%,'invert')
            axis(ax,'image')
            %% Assignment
            in.x = X;
            in.z = Z;
            in.mask = mask;
            assignin('base','velocity',in)
            
            %% Functions
            function viewingMode(type,u,v,im)
                %% DISCARDED NESTED FUNCTION - BUT HAS EXAMPLES OF OTHER FIELDS TO PLOT
                %% Calculations
                
                % calculate the velocity magnitude
                u_mag=(u.^2+v.^2).^0.5;
%                 u_max=max(max(u_mag));
%                 u_mag=u_mag/u_max;
                
                % calculate vorticity
                vor=vorticity(u, v);
                vor_max=max(max(abs(vor)));
                vor=vor/vor_max;
                
                
                % calculate the 2nd invariant
                Q=invariant2_factor(u, v, 1, 1);
                %% Plotting
                [x,y]=meshgrid(X,Z);
                [sx,sy]=meshgrid(imresize(X,1/10),imresize(Z,1/10));
                switch type
                    case 'none'               
                        imagesc(X,Z,im)
                        hold on
                        colormap gray
                        
                        h=streamslice(X,Z,u,v,2);
                        set(h, 'Color', 'r');
                        hold off
                        
                    case 'vectors'
                        imagesc(X,Z,im)
                        hold on
                        colormap gray
                        h = vis_flow(u,v,...
                            40, 1, 3, 'm',ds);
                        
                        set(h, 'Color', 'r');
                        hold off
                        
                    case 'velocity'
                        % plot velocity magnitude field
%                         ulims=[0, 1];
                        imagesc(X,Z,u_mag)%,ulims);
                        xlabel('x (pixels)');
                        ylabel('y (pixels)');
                        axis normal xy;
                        
                        title('Velocity Magnitude Field');
                        colorbar;
                        hold on;
                        
                        % plot streamlines
                        
                        h=streamslice(x,y,u,v,2);
                        set(h, 'Color', 'yellow');
                        hold off;             
                    case 'vorticity'
                        % plot Vorticity field
                        vlims=[-1, 1];
                        imagesc(X,Z,vor,vlims);
                        xlabel('x (pixels)');
                        ylabel('y (pixels)');
                        
                        
                        title('Vorticity Field');
                        colorbar;
                        colormap jet;
                        hold on;
                        
                        % plot streamlines
                        h=streamslice(x, y, u, v, 4);
                        set(h, 'Color', 'black');
                        hold off;
                        axis normal xy;
                    case 'Q'
                        %%
                        Qlims=[0, 0.1];
                        imagesc(Q,Qlims);
                        xlabel('x (pixels)');
                        ylabel('y (pixels)');
                        axis image;
                        set(gca,'YDir','reverse');
                        title('Q Field');
                        colorbar;
                end
            end
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'ax'
                            ax = namevalue;
                        case 'input'
                            in = namevalue;
                        case {'Image','u','v','U','Vorticity'}
                            view = items{m};
                        case {'None','Vectors','Streamlines'}
                            overlay = items{m};
                        case 'mask'
                            mask = true;
                    end
                    m = m+1;
                    if m>length(items);break;end
                end
                function out = namevalue
                    out = items{m+1};
                    m   = m+1;
                end
            end
        end
        function out    = getVelocity(obj,imageindex,varargin)
            %% GETVELOCITY is a function that calculates Optical Flow velocity fields from data
            %_________________________________________________________________
            %_________________________________________________________________
            % Input arguments:
            %   REQUIRED:
            %   imageindex  -   [int] 
            %                   Index of the time-series to be loaded.
            %   
            %   OPTIONAL:
            %   Argument parsing: (argval)
            %       dt      -   [int] 
            %                   Interval between image indexes for Optical
            %                   Flow algorithm (imageindex, imageindex+dt)
            %                   DEFAULT IS 1.
            %
            %       method  -   [string] 
            %
            %                   'fast'; Horn-Schunck method (bounded
            %                   iterative scheme)
            %
            %                   'slow'; Lie-Shun method (residual
            %                   minimizing scheme, uses Horn-Schunck
            %                   method, roughly 10-20x slower than 'fast'
            %                   method).
            %
            %       'mask'  -   [string]
            %                   Masking state;
            %                   true OR 1 to apply mask [DEFAULT]
            %                   false OR 0 otherwise
            %
            %   Argument parsing: (argname,argval)
            %   ARGNAME     |   ARGVAL
            %   ------------|-------------
            %   ~~~~~~~~~~~~~~~~~~~~~~~ GENERIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %   'resize'    -   Image resizing factor for variable
            %                   performance. Accepts values between
            %                   (0,inf). 
            %                   DEFAULT IS 1 (i.e., no resizing).
            %
            %   'alpha'     -   Convolution kernel size. 
            %                   DEFAULT IS 20.
            %   
            %   'loop'      -   [1xn, int] (UNCODED)
            %                   Vector for sequencing a for-loop.
            %                   DEFAULT IS EMPTY.
            %
            %   'mask'      -   [bool]
            %                   Masking state;
            %                   true OR 1 to apply mask [DEFAULT]
            %                   false OR 0 otherwise
            %   
            %   ~~~~~~~~~~~~~~~~~~~ HORN-SCHUNCK ("fast") ~~~~~~~~~~~~~~~~~~~
            %   'iterations'-   Iteration limit.
            %                   DEFAULT IS 200.
            %
            %   ~~~~~~~~~~~~~~~~~~~~~ LIE-SHUN ("fast") ~~~~~~~~~~~~~~~~~~~~~
            %   *(I'm not precisely sure what all of these do)
            %   'lambda'    -   "Regularization parameter".
            %                   DEFAULT IS 100.
            %   
            %   'illumination'- Local illumination intensity adjustment.
            %                   DEFAULT IS 0 (bypass).
            %      
            %   'filter'    -   Gaussian filter size for removing random
            %                   noise in images.
            %                   DEFAULT IS 5.   
            %_________________________________________________________________
            %_________________________________________________________________
            %   Example:
            %   % Preamble
            %   data        = exptdata(filename);
            %   imageindex  = 100;
            %   fig         = figure;
            %   ax          = axes(fig);
            %
            %   out = data.getVelocity(imageindex,'method','slow','mask');
            %_________________________________________________________________
            %_________________________________________________________________
            %% Inputs
            dt          = 1;
            method      = 'fast';
            alpha       = 20;
            lambda      = 100;
            illumination= 0;
            filter      = 5;
            iterations  = 2e2;
            resize      = 1;
            timeseries  = false;
            mask        = false;
            loop        = [];
            parseInput(varargin)
            params      = struct;

            if isempty(loop)
                loop = imageindex;
            end
            %% Implementation
            j = 1;
            for imageindex=loop
                %% Parsing
                im1     = obj.loadImage(imageindex);
                im2     = obj.loadImage(imageindex+dt);
                if resize~=1
                    im1 = imresize(im1,resize);
                    im2 = imresize(im2,resize);
                end

                if timeseries&&imageindex==loop(1)
                    % only interested in strfcn right now
%                     strfcn = zeros(size(im1,1),numel(loop));
                end

                %% Methods
                % undo mask
                im1(isnan(im1))=0;
                im2(isnan(im1))=0;
                out     = velocity_field(gather(im1),gather(im2),method,alpha,iterations,lambda,illumination,filter);
                

                %% Output
                ds      = obj.dimensions.ds;
                vthres  = 0e-5;
                convert = obj.dimensions.dt/ds;
                out.u  = out.u*convert;
                out.u  = out.u.*(abs(out.u)>vthres);
                out.v  = out.v*convert;
                out.v  = out.v.*(abs(out.v)>vthres);
                out.idx = imageindex;
                out.didx= dt;

                % masking
                if mask
                    out.u = obj.applyMask(u,imageindex);
                    out.v = obj.applyMask(v,imageindex);
                end

                obj.velocity = out;
                if timeseries
                    %                     strfcn(:,j) = obj.getStreamfunction(imageindex);
                    %                     j = j+1;
                    %                     displayProgress('Processing',imageindex,loop(1),loop(end));
                end
            end
            if timeseries
                warning('Time series functionality not yet coded')
%                 out.strfcn = strfcn;
%                 out.t = obj.times(loop)-obj.times(1);
            end
            params.method       = method;
            params.alpha        = alpha;
            params.iterations   = iterations;
            params.resize       = resize;
            out.params          = params;
            obj.velocity = out;
            %% Functions
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'lambda'
                            lambda = namevalue;
                        case 'illumination'
                            illumination = namevalue;
                        case 'filter'
                            filter = namevalue;
                        case 'alpha'
                            alpha = namevalue;
                        case 'iterations'
                            iterations = namevalue;
                        case {'slow','fast'}
                            method = items{m};
                        case {'resize','rescale'}
                            resize = namevalue;
                        case 'mask'
                            mask    = true;
                        case {'ts','timeseries'}
                            loop = namevalue;
                            if ~isempty(loop)
                                timeseries = true;
                            end
                        otherwise
                            if isa(items{m},'double')
                                dt = round(items{m});
                            end
                    end
                    m = m+1;
                    if m>length(items);break;end
                end
                function out = namevalue
                    out = items{m+1};
                    m   = m+1;
                end
            end
        end
        function out    = applyMask(obj,in,imageindex)
            %% APPLYMASK is a function that masks velocity fields with NaN data from the original image
            %_________________________________________________________________
            %_________________________________________________________________
            % Input argument:
            %   REQUIRED:
            %   in          -   [nxm matrix, double]
            %                   Input velocity field matrix.
            %
            %   OPTIONAL:
            %   imageindex  -   [int] 
            %                   index of the time-series to be loaded
            %
            %
            %_________________________________________________________________
            %_________________________________________________________________
            % Output argument:
            %   out         -   [nxm matrix, double] 
            %                   Masked input velocity field matrix.
            %_________________________________________________________________
            %_________________________________________________________________
            %%
            if nargin<3
                imageindex = 1;
            end

            mask = ~isnan(obj.loadImage(imageindex));
            mask = imresize(mask,size(in));
            in(~mask)=0;
            out = in;
            %                 noslipmask = imgaussfilt(double(mask),5);
            %                 ux = ux.*noslipmask;
            %                 uy = uy.*noslipmask;

        end
    end
    methods % IMAGE LOADING AND DIMENSIONS
        function im = loadImage(obj,imageindex)
            % LOADIMAGE loads the 'SST' with [1 1 idx] x-z-t index,[Inf Inf 1] x-z-t range
            %_________________________________________________________________
            %_________________________________________________________________
            %   Input argument:
            %   imageindex  -   [int] 
            %                   index of the time-series to be loaded 
            %_________________________________________________________________
            %_________________________________________________________________
            %   Output argument:
            %   im          -   [nxm matrix, double]
            %                   Image data.
            %_________________________________________________________________
            %_________________________________________________________________
            %% Reading
            % Accept appropriate indexes
            if imageindex<=0
                imageindex = 1;
            elseif imageindex>numel(obj.times)
                imageindex = numel(obj.times);
            end
            im = ncread(obj.ncinfo.Filename,'SST',[1 1 imageindex],[Inf Inf 1]);
            %% Manipulation
            im(im==0)=NaN; % Set 0's to NaN
            %             [dx,dy] = gradient(im);
            %             im = dy.^2+dx.^2;
            obj.image = im; % Save current image output to global properties
        end
        function [x,z,t,ds,dt] = getDimensions(obj,centered)
            %% GETDIMENSIONS is a function that returns x,z,t vectors and ds,dt values
            %_________________________________________________________________
            %_________________________________________________________________
            %   Input argument:
            %   centered    -   [bool] 
            %                   Set [x,z] vectors to be symmetric about 0.
            %                   DEFAULT IS FALSE.
            %_________________________________________________________________
            %_________________________________________________________________
            %_________________________________________________________________
            % Output arguments:
            %   'x'         -   [1xn vec, double] 
            %                   X-wise vector in meters.
            %
            %   'z'         -   [1xn vec, double] 
            %                   Z-wise vector in meters.
            %
            %   't'         -   [1xn vec, duration] 
            %                   T-wise vector in duration format.
            %                   (Use "seconds(t)" or similar to return to
            %                   double) 
            %
            %   'ds'        -   [double]
            %                   Pixel to meter ratio. 
            %
            %   'dt'        -   [duration]
            %                   Image index time difference.
            %_________________________________________________________________
            %_________________________________________________________________
            %   Example:
            %   % Premable:
            %   data = exptdata(filename);
            %   
            %   [x,z,t,ds,dt] = data.getDimensions;
            %% Conversion and output
            ds      = obj.dimensions.ds;
            if isempty(obj.image)
                obj.loadImage(1);
            end
            im      = obj.image;
            z       = [0:size(im,1)-1]/ds; %#ok<*NBRAK>
            x       = [0:size(im,2)-1]/ds;
            t       = obj.times-obj.times(1);
            dt      = mean(diff(t));

            if nargin==2
                if centered
                    x = x-max(x)/2;
                    z = z-max(z)/2;
                end
            end
        end
    end
end