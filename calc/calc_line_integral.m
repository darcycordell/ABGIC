function [gic1d,gic3d] = calc_line_integral(lines,tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,method)
%Calculate line integral along transmission line paths given a start point
%and end point of the transmission lines
%
% Inputs: lines is a cell of matrices. Each cell entry represents one
%               transmission line and each tranmission line consists of a
%               (2 x N) matrix of latitude and longitude pairs. The lines should be
%               pre-processed so that the distance between points on the
%               tranmission line is small enough to be valid over the
%               interpolation (usually N ~ 100).
%         tind is the time index (or indices) where you want to compute the
%               line integral
%         d is a data structure containing the locations of MT impedances
%         ex3d, ey3d are (nt by ns) matrices of electric fields (V/km) where nt is
%               the number of time steps and ns is the number of MT sites
%         ex1d, ey1d are (nt by ngrid) matrices of electric field (V/km) where nt
%               is the number of time steps and ngrid is the number of grid
%               points where the e-field was interpolated
%         LAT and LON are nx by ny meshgrids where the 1D efield was
%                interpolated
%         method is the interpolation method (e.g. cubic, nearest neighbor)
%
%
% Outputs:
%       gic1d is the (length(tind) x ntran) matrix of gic calculations
%               using the electric field calculated using 1-D impedances
%               one calculation for each time index and each transmission line
%       gic3d is the (length(tind) x ntran) matrix of gic calculations
%               using the electric field calculated using 3-D MT impedances

%%

warning('off');

%Initialize large matrices
gic1d = zeros(length(tind),length(lines)); gic3d = zeros(size(gic1d));

%The origin of the survey used for lat-long to UTM conversions
origlon =(max(d.loc(:,2))-min(d.loc(:,2)))/2+min(d.loc(:,2));
origlat = (max(d.loc(:,1))-min(d.loc(:,1)))/2+min(d.loc(:,1));
nlines = size(lines,2);

for sidx = 1:nlines %Loop through lines
    
    g3d = zeros(length(tind),1);
    g1d = zeros(length(tind),1);
    
    %Line vertices in longitute (rx) and latitude (ry)
    rx = lines{sidx}(1,:);
    ry = lines{sidx}(2,:);
    %Convert lat long to metres (UTM). Note: y = easting; x = northing 
    % Two methods: Tested both and it doesn't change the results at all.
    %[y,x] = geo2utm(ry,rx,ry(round(length(ry)/2)),rx(round(length(rx)/2))); convert using the transmission line segment midpoint as a reference
    [y,x] = geo2utm(ry,rx,origlon,origlat); %convert using the full dataset midpoint as an origin
    tic
    
    %Interpolate transmission line values and convert to V/m
    % Here I use scatteredInterpolant to generate 4 functions for ex1d,
    % ey1d, ex3d, and ey3d. Once I have the interpolant function, then I
    % enter the loop and just replace the function values since d.loc never
    % changes. This is faster than e.g. gridding on every time step.

    Fex3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ex3d(tind(1),:).',method);
    exnn3d = Fex3d(rx,ry);
    
    Fey3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ey3d(tind(1),:).',method);
    eynn3d = Fey3d(rx,ry);
    
    g3d(1) = abs(nansum(abs(diff(x)).*(exnn3d(1:end-1)+exnn3d(2:end))/2+abs(diff(y)).*(eynn3d(1:end-1)+eynn3d(2:end))/2));

    Fex1d = scatteredInterpolant(LAT(:),LON(:),ex1d(tind(1),:).',method);
    exnn1d = Fex1d(rx,ry);
    
    Fey1d = scatteredInterpolant(LAT(:),LON(:),ey1d(tind(1),:).',method);
    eynn1d = Fey1d(rx,ry);
    
    g1d(1) = abs(nansum(abs(diff(x)).*(exnn1d(1:end-1)+exnn1d(2:end))/2+abs(diff(y)).*(eynn1d(1:end-1)+eynn1d(2:end))/2));
    
    for tidx = 2:length(tind) %Loop through each time step

        Fex3d.Values = ex3d(tind(tidx),:).';
        exnn3d = Fex3d(rx,ry);
        
        Fey3d.Values = ey3d(tind(tidx),:).';
        eynn3d = Fey3d(rx,ry);
        
        Fex1d.Values = ex1d(tind(tidx),:).';
        exnn1d = Fex1d(rx,ry);
        
        Fey1d.Values = ey1d(tind(tidx),:).';
        eynn1d = Fey1d(rx,ry);
        
        
        %Perform line integral (x and y are in meters, e is in V/m)
        %Take the absolute value because we do not care so much about
        %whether GIC is positive or negative but only the absolute
        %difference between methods
        % Resources:
        %       https://www.mathworks.com/matlabcentral/answers/441416-numerical-calculation-of-line-integral-over-a-vector-field
        %       https://ocw.mit.edu/ans7870/18/18.013a/textbook/HTML/chapter25/section04.html
        g3d(tidx) = (nansum(abs(diff(x)).*(exnn3d(1:end-1)+exnn3d(2:end))/2+abs(diff(y)).*(eynn3d(1:end-1)+eynn3d(2:end))/2));

        g1d(tidx) = (nansum(abs(diff(x)).*(exnn1d(1:end-1)+exnn1d(2:end))/2+abs(diff(y)).*(eynn1d(1:end-1)+eynn1d(2:end))/2));
        
        if rem(tidx,3600)==0
            disp(['Hour #',num2str(tidx/3600),' Complete'])
        end

    end
    
    gic3d(:,sidx) = g3d;
    gic1d(:,sidx) = g1d;
toc
    disp(['Transmission Line #',num2str(sidx),' of ',num2str(nlines),' Completed...........................'])

end


warning('on')


%%
% For debugging to make sure the line integral is working correctly

% figure(30);
% u.dx = 0.1; u.dy = 0.1;
% xgrid = d.lim(1):u.dx:d.lim(2);     
% ygrid = d.lim(3):u.dy:d.lim(4); 
% [X,Y] = meshgrid(xgrid,ygrid);
% 
% rhogrid = griddata(d.loc(:,2),d.loc(:,1), squeeze(ey3d(tind,:)),X,Y,method);
% pcolor(X,Y,rhogrid); hold on
% plot(d.loc(:,2),d.loc(:,1),'.k')
% 
% plot(ry,rx,'-k');
% 
% quiver(d.loc(:,2),d.loc(:,1),ey3d(tind,:)',ex3d(tind,:)','-k')
% quiver(LON(:),LAT(:),ey1d(tind,:)',ex1d(tind,:)','-r')
% 
% quiver(ry,rx,eynn1d*1000,exnn1d*1000,'-b','AutoScale','off')
% quiver(ry,rx,eynn3d*1000,exnn3d*1000,'-g','AutoScale','off')
% 
% colorbar
% caxis([-0.2 0.2])
% axis([min(ry)-0.1 max(ry)+0.1 min(rx)-0.1 max(rx)+0.1])


        
    