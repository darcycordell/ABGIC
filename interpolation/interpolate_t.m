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
% Option to set lim = NaN and interpolate on the (dlat,dlon) points instead
%

%Number of time steps in the data files
N = length(a(1).x);

%Get data matrices for interpolation
s_lat = [a(:).lat];
s_lon = [a(:).lon];
x = [a(:).x];
y = [a(:).y];


%% interpolate onto mesh

if ~isnan(lim)

    [lon,lat] = meshgrid(lim(3):dlon:lim(4),lim(1):dlat:lim(2));

else
    
    lon = dlon;
    lat = dlat;
end


t_ind = 1;
xint = NaN(size(lon,1),size(lon,2),N); % pre allocate interpolated fields
yint = NaN(size(lon,1),size(lon,2),N);

Fx_int = scatteredInterpolant(s_lon',s_lat',x(t_ind,:)','natural');
Fy_int = scatteredInterpolant(s_lon',s_lat',y(t_ind,:)','natural');

tic
Ax=func2mat(@(D) func_func2mat(D,Fx_int,lon,lat), x(t_ind,:)' );
xint=Ax*x';

Ay=func2mat(@(D) func_func2mat(D,Fy_int,lon,lat), y(t_ind,:)' );
yint=Ay*y';
toc

xint = reshape(xint,size(lat,1),size(lat,2),N);
yint = reshape(yint,size(lat,1),size(lat,2),N);

% 
% for j = 1:1:N
%     
%     Fx_int.Values = x(t_ind,:)'; 
%     Vxq = Fx_int(lon,lat);
%     xint(:,:,t_ind) = Vxq;
%    
%     Fy_int.Values = y(t_ind,:)';
%     Vyq = Fy_int(lon,lat);
%     yint(:,:,t_ind) = Vyq;
%     
%     if rem(t_ind,3600)==1
%         disp(['Time Step: ',num2str(t_ind),'; ',datestr(a(1).times(t_ind))])
%     end
%        
%     t_ind = t_ind+1;
%     j
% end
% toc
