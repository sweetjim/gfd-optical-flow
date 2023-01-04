function play(M,u,v,umean,vmean)
set(gca,'Color','none')
for i=1:size(M,3)
    imagesc(M(:,:,i),'AlphaData',~isnan(M(:,:,i)))
    if i==1
        clim=[20 30];%caxis;
    end
    if nargin<4
        hold on
        try
            vis_flow(u(:,:,i),v(:,:,i))
        catch
            vis_flow(u,v)
        end
        hold off
    elseif nargin>3
        hold on
        vis_flow(u(:,:,i),v(:,:,i),'c','r')
        vis_flow(umean,vmean,'c','w')
        hold off
    end
    caxis(clim);
    try
        cmocean('balance','pivot',0);
    end
    pause(.1)
end


