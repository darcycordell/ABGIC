function [zn, Clat, Clon, txt, nb] = load_trich_zones

    %These Zones were given by Larisa Trichtchenko and Dave Boteler from their
    %2019 report (see Figure 4.1.1 and 4.2.1 in their report).
    ab_zones = jsondecode(fileread('Alberta.geojson'));
    bc_zones = jsondecode(fileread('BritishColumbia.geojson'));
    
    ab_zone_key = [7 15 6 5 10 8 3 1 9 2 4];
    bc_zone_key = [16 17 18 19 12 20 21 22 23 13 14 11 24];

    %Load AB Zones
    zn = []; Clon = []; Clat = [];
    for i = 1:length(ab_zones.features)
        lines = squeeze(ab_zones.features(i).geometry.coordinates);
        [zn(ab_zone_key(i)).lon,zn(ab_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
    end

    nb = length(zn);
    %Load BC Zones
    for i = 1:length(bc_zones.features)
        lines = squeeze(bc_zones.features(i).geometry.coordinates);
        [zn(bc_zone_key(i)).lon,zn(bc_zone_key(i)).lat] = poly2cw(lines(:,1),lines(:,2));
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