function cmap = cmapCheck(cmap)
switch colormapclasses(cmap)
    case 'cmocean'
        cmap = cmocean(cmap);
    otherwise
        cmap = colormap(cmap);
        return
end
end