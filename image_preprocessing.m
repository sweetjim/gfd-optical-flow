function output = image_preprocessing(input,type,varargin)
output = input; % default is no pre-processing
if nargin==1
    return
end
resolution = max(size(input));
parseInput(varargin)

switch type
    case 'cartesian'
        return
    case {'polar','ridge','square'}
        %% Convert to polar coordinates
        [r,c] = size(input);
        clip = (diff([r c])/2);

        if clip~=0
            if clip<0
                % clipping
                zlims = -clip:r+clip-1;
                xlims = 1:c;
            elseif clip>0
                % clipping
                zlims = 1:r;
                xlims = clip:r-clip-1;
            end
            % padding
            input = padarray(input,.5*abs([r,c]-max([r,c])));
            input(input==0)=nan;

            % if clipping
            %input = input(zlims,xlims);
            if strcmp(type,'square')
                output = input;
                return
            end
        end
        dtheta  = resolution;
        dr      = resolution;

        [output,r,theta] = ImToPolar(input,...
            .36,...         % radius inner start [0,1]
            .8,...          % radius outer end   [0,1]
            dr,dtheta);
        if strcmp(type,'polar')
            return
        end
        %% Set to ridge following coordinates
    case 'polarinv'
        output = PolarToIm(input, ...
            .36, ...
            .8, ...
            size(output,1), ...
            size(output,2));
end
    function parseInput(varargin)
    m = 1;
    items = varargin{:};
    for k=1:length(items)
        switch items{m}
            case {'res','resolution'}
                val = namevalue;
                if isempty(val)
                    val = resolution;
                end
                resolution = val;
            case 'coarse'
                resolution = 100;
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

