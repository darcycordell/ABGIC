function [d,in,indzones] = load_assign_impedance(zn,mtfile)
%
% Function which loads MT impedance data and assigns each MT site to the
% appropriate 1-D zone from Trichtchenko et al. (2019).
%
% Inputs:
%       zn: The zones (from Trichtchenko et al., 2019) are loaded into the "zn"
%           structure containing lat/long points for the polygons for each zone
%
%       mtfile: a mat file containing the MT data in a "d" structure
%               mtfile must be on path or in current directory
%
% Outputs:
%       d: a data structure containing all impedance info
%               Note that loaded impedances are in SI units (Ohm)
%               defined by E = ZH (E = V/m and H = A/m)
%       in: a logical mask which assigns each MT site to a specific zone
%               size in = (# of MT sites) x (# of zones)
%       indzones: the index of the first site in each zone
%

[zn, Clat, Clon, txt, nb] = load_trich_zones;
val = load(mtfile);
d = val.d;

% plot_lim = [47.5 61 -125.5 -109];
% provinces = shaperead('province.shp','UseGeoCoords',true);
% states = shaperead('usastatehi','UseGeoCoords',true);
% initialize_map(plot_lim,zn,provinces,states,101,true)
% textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')


%Find which sites are in each zone polygon
nb = length(zn);
in = false(d.ns,nb);
for i = 1:nb
    in(:,i) = inpolygon(d.loc(:,1),d.loc(:,2),zn(i).lat,zn(i).lon);
    %geoshow(d.loc(in(:,i),1),d.loc(in(:,i),2),'displaytype','point','marker','o','markerfacecolor','red','markeredgecolor','black','markersize',4)
end

% geoshow(Clat,Clon,'displaytype','point','marker','o','markerfacecolor','yellow','markeredgecolor','black','markersize',4)
% textm(Clat,Clon,txt,'HorizontalAlignment','center')

%There are some MT sites in the survey which are NOT inside any polygons
%(e.g. a couple American sites and a couple Sask sites)
ind = find(all(in==0,2)); %This finds the site indices which are not in a polygon

%Loop through the sites that aren't in a polygon and assign them to a
%polygon based on the nearest polygon
for i = 1:length(ind)
    %%
    for k = 1:nb
        %Find minimum distance from site outside polygon to edge of all
        %polygons
        [mn(k)] = min(distance(d.loc(ind(i),1),d.loc(ind(i),2),zn(k).lat,zn(k).lon));
    end
    
    %Find index of the polygon which has an edge closest to the site
    [~,id] = min(mn);
    
    %Set the logical mask as true at that point instead of false.
    in(ind(i),id) = true;
    
    %For debugging purposes to check sites are being assigned correctly
     %geoshow(d.loc(ind(i),1),d.loc(ind(i),2),'displaytype','point','marker','o','markerfacecolor','red','markeredgecolor','black','markersize',4)
     %geoshow(Clat(id),Clon(id),'displaytype','point','marker','o','markerfacecolor','yellow','markeredgecolor','black','markersize',4) 
     %plotm([Clat(id) d.loc(ind(i),1)],[Clon(id) d.loc(ind(i),2)],'-y','LineWidth',2) 
end

%There are some points for which this algorithm doesn't work. I could've
%made it more sophisticated by using e.g. line segment calculations and
%line intersections, but I just left it as is. There are 7 sites which are
%assigned to the wrong zone. Here I assign the logical mask manually:
in(ind(4),1) = false; in(ind(4),13) = true; %id = 21;
            %in(ind(12),19) = false; in(ind(12),23) = true; %id = 23;
            %in(ind(15),16) = false; in(ind(15),23) = true; %id = 23;
in(ind(19),3) = false; in(ind(19),4) = true; %id = 11;
in(ind(20),3) = false; in(ind(20),4) = true; %id = 11;
in(ind(21),3) = false; in(ind(21),4) = true; %id = 11;
in(ind(22),3) = false; in(ind(22),4) = true; %id = 11;


%Site NAB875 is in the polygon labelled "Zone 2" but in the original
%Trichtchenko report, this polygon is not given a zone # or resistivity
%model.
in(450,15) = false; in(450,7) = true; %513 Site Cedar data has NAB875 at index 451

%%
%For debugging purposes, you can uncomment and check each zone to ensure
%that points are in the right spots
% zo = 24;
% plot_lim = [47.5 61 -125.5 -109];
% provinces = shaperead('province.shp','UseGeoCoords',true);
% states = shaperead('usastatehi','UseGeoCoords',true);
% initialize_map(plot_lim,zn,provinces,states,101,true)
% textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')
% geoshow(d.loc(in(:,zo),1),d.loc(in(:,zo),2),'displaytype','point','marker','o','markerfacecolor','yellow','markeredgecolor','black','markersize',4)
%%
indzones = find(~all(in==0));

%For debugging purposes
% for i = indzones
%     geoshow(d.loc(in(:,i),1),d.loc(in(:,i),2),'displaytype','point','marker','o','markerfacecolor',rand(3,1),'markeredgecolor','black','markersize',4)
% end