function previewExpt(callingApp,resolution)
if nargin==1
    resolution = [];
end
OE          = callingApp.openExpt;
if iscell(OE)
OE = OE{1};
end
file2read   = OE;
% file2read   = fullfile(OE.folder,OE.name);
idx         = callingApp.imIndex;
%%
try
    im = ncread(file2read,'SST',[1 1 idx],[Inf Inf 1]);
catch
    im = ncread(file2read,'u',[1 1 idx],[Inf Inf 1]);
end
%% image manipulation
im(im==0)=NaN;
im=image_preprocessing(im,'polar','res',resolution);
%%
caxis([20 25])
imagesc(im,'AlphaData',~isnan(im))
caxis([23 27])
set(gca,'Colormap',1-cmocean('thermal'))
% colorbar
assignin("base","im",im)
end

