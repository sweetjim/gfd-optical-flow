function addlabels(varargin)
%% Add labels to an axes


%%
ax          = [];
x_str       = '';
y_str       = '';
z_str       = '';
title_str   = '';
fs          = 10;
latex       = false;
parseInput(varargin);

if isempty(ax)
    ax = gca;
end

for i=1:numel(ax)
    axi = ax(i);
    if numel(title_str)>1 && isa(title_str,'cell')
       title_stri = title_str{i};
    else
        title_stri = title_str;
    end
    if latex
        try
            xlabel(axi,x_str,'Interpreter','latex')
            ylabel(axi,y_str,'Interpreter','latex')
            title(axi,title_stri,'Interpreter','latex')
            set(axi,'TickLabelInterpreter','latex')
            zlabel(axi,z_str,'Interpreter','latex')
        end
    else
        try
            xlabel(axi,x_str)
            ylabel(axi,y_str)
            title(axi,title_stri)
            zlabel(axi,z_str)
        end
    end
end
try
    set(ax,'FontSize',fs)
catch
    if isa(ax,'matlab.graphics.layout.TiledChartLayout')
        tiletitle = get(ax,'Title'); set(tiletitle,'FontSize',fs);
        tilex = get(ax,'Xlabel'); set(tilex,'FontSize',fs);
        tiley = get(ax,'Ylabel'); set(tiley,'FontSize',fs);
    end
end


%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'latex'
                    latex   = true;
                    %% Name-value arguments
                case 'ax'
                    ax      = namevalue;
                case 'title'
                    title_str   = namevalue;
                case {'x','x_str'}
                    x_str   = namevalue;
                case {'y','y_str'}
                    y_str   = namevalue;
                case {'z','z_str'}
                    z_str   = namevalue;
                case 'fs'
                    fs      = namevalue;
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