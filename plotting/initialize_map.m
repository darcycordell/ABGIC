function initialize_map(plot_lim,zn,provinces,states,fignum,plot_zone)
%%
screensize=get(groot,'Screensize');
fig=figure(fignum);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
%%
if plot_zone
    %Load and Plot AB Zones
    for i = 1:length(zn)
        geoshow(zn(i).lat,zn(i).lon,'displaytype','polygon','facealpha',0.5); hold on
    end
end

geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0)
geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0)