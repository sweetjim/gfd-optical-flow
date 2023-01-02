function previewExpt(callingApp)
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

[r,c] = size(im);
clip = (diff([r c])/2);
if clip~=0
    if clip<0
        zlims = -clip:r+clip-1;
        xlims = 1:c;
    elseif clip>0
        zlims = 1:r;
        xlims = clip:r-clip-1;
    end
    im = im(zlims,xlims);
end
%% image manipulation
im(im==0)=NaN;
%%
caxis([20 25])
% imagesc(ImToPolar(im,0.37,1,size(im,1),size(im,1)))
imagesc(im,'AlphaData',~isnan(im))
caxis([23 27])
set(gca,'Colormap',1-cmocean('thermal'))
% colorbar
assignin("base","im",im)
end

