function [xint,yint] = interpolate_e_from_Z(d,x,y,lat,lon)
% Loads MT impedance data locations and ex and ey time domain E-field data
% at those locations and then interpolates onto a regular grid defined by
% LAT, LON
%
%   Inputs:
%       d = data structure from MT with ns sites
%       x, y = matrices of time-domain E-field data with size (nt x ns)
%               where nt is the number of time steps
%       lat, lon = latitude and longitude coordinates defined on a meshgrid
%                   with size (nx x ny). Total number of grid points is
%                   N = nx*ny
% 	Output matfile includes:
%       xint,yint: Interpolated E fields in V/m with size (nt x N)
%

%Number of time steps in the data files
N = size(x,1);

%Get data matrices for interpolation
s_lat = d.loc(:,1);
s_lon = d.loc(:,2);

%% interpolate onto mesh

t_ind = 1;
xint = NaN(size(lat,1),size(lon,2),N); % pre allocate interpolated fields
yint = NaN(size(lat,1),size(lon,2),N);

Fx_int = scatteredInterpolant(s_lon,s_lat,x(t_ind,:)','natural');
Fy_int = scatteredInterpolant(s_lon,s_lat,y(t_ind,:)','natural');

tic
Ax=func2mat(@(D) func_func2mat(D,Fx_int,lon,lat), x(t_ind,:)' );
xint=Ax*x';

Ay=func2mat(@(D) func_func2mat(D,Fy_int,lon,lat), y(t_ind,:)' );
yint=Ay*y';
toc

%xint = reshape(xint,size(lat,1),size(lat,2),N);
%yint = reshape(yint,size(lat,1),size(lat,2),N);

xint = xint';
yint = yint';

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
