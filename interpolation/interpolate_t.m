function [xint,yint,lon,lat] = interpolate_t(a,lim,dlon,dlat)
% Load B field data from multiple stations
% Interpolate B field data over time and space
% Save interpolated B field into mat file
%
%   Inputs:
%       a = data structure from load_mag_data
%       lim = [minlat, maxlat, minlon, maxlon]
%       dlon, dlat = cell sizes in decimal degrees
% 	Output matfile includes:
%       xint,yint: Interpolated B fields in nT (nlat x nlon x ntimesteps)
%       lat: meshgrid of latitudes (nlat x nlon)
%       lon: mesghrid of longitudes (nlat x nlon)
%
% **Note that mag data files must cover the exact same time range**
%
%

%Number of time steps in the data files
N = length(a(1).x);

%Build data matrices for interpolation
s_lat = [];s_lon = [];x=[];y=[];%z=[];
for i = 1:length(a)
    s_lat = [s_lat a(i).lat];
    s_lon = [s_lon a(i).lon];
    x = [x a(i).x];
    y = [y a(i).y];
    %no z interpolation
end


%% interpolate onto mesh

[lon,lat] = meshgrid(lim(3):dlon:lim(4),lim(1):dlat:lim(2));

t_ind = 1;
xint = NaN(size(lon,1),size(lon,2),N); % pre allocate interpolated fields
yint = NaN(size(lon,1),size(lon,2),N);
tic
for j = 1:1:N
    
    Fx_int = scatteredInterpolant(s_lon',s_lat',x(t_ind,:)');
    Vxq = Fx_int(lon,lat);
    xint(:,:,t_ind) = Vxq;

    Fy_int = scatteredInterpolant(s_lon',s_lat',y(t_ind,:)');
    Vyq = Fy_int(lon,lat);
    yint(:,:,t_ind) = Vyq;
    
    if rem(j,3600)==1
        disp(['Time Step: ',num2str(j),'; ',datestr(a(1).times(j))])
    end
       
    t_ind = t_ind+1;
end
toc