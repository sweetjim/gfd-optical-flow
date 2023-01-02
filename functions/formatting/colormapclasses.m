function maps = colormapclasses(varargin)
cmoceanmaps = {'thermal','balance','haline','delta','solar','curl'...
    'ice','diff','gray','tarn',...
    'oxy','deep','dense','phase',...
    'algae','matter','turbid','topo',...
    'speed','amp','tempo','rain'};
builtinmaps = {'parula','jet','hsv','hot','cool','spring',...
    'spring','summer','autumn','winter',...
    'gray','bone','copper','pink','lines',...
    'colorcube','prism','flag','white'};

if nargin>0
    map = varargin{1};
    maps = 'builtin';
    if contains(map,cmoceanmaps)
        maps = 'cmocean';
    end
    return
end

maps = cat(2,cmoceanmaps,builtinmaps);

end

