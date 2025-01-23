function [zn, Clat, Clon, txt, nb] = load_trich_zones

    %These Zones were given by Larisa Trichtchenko and Dave Boteler from their
    %2019 report (see Figure 4.1.1 and 4.2.1 in their report).
    ab_zones = jsondecode(fileread('Alberta.geojson'));
    bc_zones = jsondecode(fileread('BritishColumbia.geojson'));
    sk_zones = jsondecode(fileread('Saskatoon.geojson'));
    mb_zones = jsondecode(fileread('Manitoba.geojson'));
    
    ab_zone_key = [7 15 6 5 10 8 3 1 9 2 4];
    bc_zone_key = [16 17 18 19 12 20 21 22 23 13 14 11 24];
    sk_zone_key = 25:40;
    mb_zone_key = 41:53;

    %Load AB Zones
    zn = []; Clon = []; Clat = [];
    for i = 1:length(ab_zones.features)
        lines = squeeze(ab_zones.features(i).geometry.coordinates);
        [zn(ab_zone_key(i)).lon,zn(ab_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
        zn(ab_zone_key(i)).province = 'AB';
    end

    nb = length(zn);
    %Load BC Zones
    for i = 1:length(bc_zones.features)
        lines = squeeze(bc_zones.features(i).geometry.coordinates);
        [zn(bc_zone_key(i)).lon,zn(bc_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
        zn(bc_zone_key(i)).province = 'BC';
    end
    nb = length(zn);

    %Load SK Zones
    for i = 1:length(sk_zones.features)
        lines = squeeze(sk_zones.features(i).geometry.coordinates);
        [zn(sk_zone_key(i)).lon,zn(sk_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
        zn(sk_zone_key(i)).province = 'SK';
    end
    nb = length(zn);

    %Load MB Zones
    for i = 1:length(mb_zones.features)
        lines = squeeze(mb_zones.features(i).geometry.coordinates);
        [zn(mb_zone_key(i)).lon,zn(mb_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
        zn(mb_zone_key(i)).province = 'MB';
    end
    nb = length(zn);

    %Find the centroid of each zone polygon using the centroid formula
    %   see: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    for i = 1:nb

        lon = zn(i).lon;
        lat = zn(i).lat;
        Sx = 0; Sy = 0; Area = 0;
        for k = 1:length(lon)-1
            sx = (lon(k)+lon(k+1))*(lon(k)*lat(k+1)-lon(k+1)*lat(k));
            Sx = Sx + sx;

            sy = (lat(k)+lat(k+1))*(lon(k)*lat(k+1)-lon(k+1)*lat(k));
            Sy = Sy + sy;

            area = lon(k)*lat(k+1)-lon(k+1)*lat(k);
            Area = Area + area;
        end

        Area = 0.5*Area;
        Clon(i) = (1/(6*Area))*Sx;
        Clat(i) = (1/(6*Area))*Sy;

    end

    for i = 1:length(Clat)
        txt{i} = ['Zone ',num2str(i)];
    end


%%

%     plot_lim = [48 61 -140 -90];
%     provinces = shaperead('province.shp','UseGeoCoords',true);
%     states = shaperead('usastatehi','UseGeoCoords',true);
%     %Figure 1a: 1-D Impedance
%     screensize=get(groot,'Screensize');
%     fig=figure(1);
%     clf
%     set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
%     worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
%     setm(gca,'FontSize',18)
%     
%     geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)
%     geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)
% 
%     cmap = hsv(nb);
%     for i = 1:nb
%         patchm(zn(i).lat,zn(i).lon,cmap(i,:));
%         textm(Clat(i),Clon(i),num2str(i))
%     end