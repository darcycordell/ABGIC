function [gic1d,gic3d,avg_strike,avg_r,std_strike,emag,simple_LV] = calc_line_integral(lines,tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,method,compute_coupling)
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


%There are two different methods to do the interpolation of the E-fields
%onto the transmission lines. The first uses a FEX script called func2mat
%which does the interpolation of *all* transmission line points and all
%time steps in a single step. This is *incredibly* fast for large problems
%and allows for the calculation of the line voltages for a 1000s or 10000s
%of points in tens of seconds.
%
%The second option does the interpolation in nested for loops, looping over
%each transmission line and each time interval. This is incredibly slow
%when considering many time steps. For example, 1000s of time steps would
%take 10s of *minutes*. However, there is a tradeoff somewhere, because
%this method is much faster if you are only doing a few timesteps. For
%example, doing only 1 timestep takes only 2 seconds, but still takes about
%35 seconds using the other method. 
%
%I tested this a bit and there doesn't seem to be any hard and fast rule
%for which is faster. But, roughly it seems that if you are doing >500 time
%steps, then the first option is faster, so I've put an if statement here
if length(tind)>=500
    speedup_option = 1;
else
    speedup_option = 0;
end


warning('off');

%Initialize large matrices
gic1d = zeros(length(tind),length(lines)); gic3d = zeros(size(gic1d));

if compute_coupling == 1 || compute_coupling == 3
    avg_strike = zeros(size(gic1d));
    avg_r = zeros(size(gic1d));
    std_strike = zeros(size(gic1d));
    emag = zeros(size(gic1d));
    simple_LV = zeros(size(gic1d));
else
    avg_strike = NaN;
    avg_r = NaN;
    std_strike = NaN;
    emag = NaN;
    simple_LV = NaN;
end

%The origin of the survey used for lat-long to UTM conversions
origlon =(max(d.loc(:,2))-min(d.loc(:,2)))/2+min(d.loc(:,2));
origlat = (max(d.loc(:,1))-min(d.loc(:,1)))/2+min(d.loc(:,1));
nlines = size(lines,2);

rx = []; ry = []; indLines = {};
for sidx = 1:nlines
    
    indLines{sidx} = 1+length(rx):length(rx)+length(lines{sidx}(1,:));
    
    rx = [rx lines{sidx}(1,:)];
    ry = [ry lines{sidx}(2,:)];
   
end

disp('Finished collecting all transmission line coordinates')

 
[y,x] = geo2utm(ry,rx,origlon,origlat); %convert using the full dataset midpoint as an origin (x = northing, y = easting)
  
Fex3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ex3d(tind(1),:).',method);
Fey3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ey3d(tind(1),:).',method);    

Fex1d = scatteredInterpolant(LAT(:),LON(:),ex1d(tind(1),:).',method);
Fey1d = scatteredInterpolant(LAT(:),LON(:),ey1d(tind(1),:).',method);



if speedup_option == 1 %"faster" option

    disp('LENGTH(TIND)>500. USING FUNC2MAT Y(:)=A*X(:) METHOD-----')
    
    disp('Interpolating E fields onto transmission line coordinates........')
    disp('Working on Ex........')
    Ax3d=func2mat(@(D) func_func2mat(D,Fex3d,rx,ry), ex3d(tind(1),:).');
    Ax1d=func2mat(@(D) func_func2mat(D,Fex1d,rx,ry), ex1d(tind(1),:).');

    disp('Finished Ex........')
    disp('Working on Ey........')

    Ay3d=func2mat(@(D) func_func2mat(D,Fey3d,rx,ry), ey3d(tind(1),:).');
    Ay1d=func2mat(@(D) func_func2mat(D,Fey1d,rx,ry), ey1d(tind(1),:).');
    
    disp('Finished Ey........')
    
    
    su_window = 3600; %to save memory, I do a loop through a set of time-steps (e.g. 1 hour)
            %I can do the matrix multiplication in one step, but this
            %results in exnn3d, eynn3d, exnn1d, and eynn1d matrices that
            %are huge and exceed 16 GB of RAM. I can save memory by doing
            %each hour at a time
    su_ind = unique([0:su_window:length(tind) length(tind)]);

    
    for su = 1:length(su_ind)-1
        
        exnn3d=Ax3d*ex3d(tind(su_ind(su)+1:su_ind(su+1)),:).';
        eynn3d=Ay3d*ey3d(tind(su_ind(su)+1:su_ind(su+1)),:).';

        exnn1d=Ax1d*ex1d(tind(su_ind(su)+1:su_ind(su+1)),:).';
        eynn1d=Ay1d*ey1d(tind(su_ind(su)+1:su_ind(su+1)),:).';
        
        for sidx = 1:nlines

            s = indLines{sidx};

            indslice = su_ind(su)+1:su_ind(su+1);

            dx = diff(x(s));
            dy = diff(y(s));

            d = distance(rx(s(1)),ry(s(1)),rx(s(end)),ry(s(end)),referenceEllipsoid('WGS84'));

            dex3 = (exnn3d(s(1:end-1),:)+exnn3d(s(2:end),:))/2;
            dey3 = (eynn3d(s(1:end-1),:)+eynn3d(s(2:end),:))/2;

            dex1 = (exnn1d(s(1:end-1),:)+exnn1d(s(2:end),:))/2;
            dey1 = (eynn1d(s(1:end-1),:)+eynn1d(s(2:end),:))/2;

            if length(s)>1
                gic3d(indslice,sidx) = nansum(repmat(dx,size(dex3,2),1).*dex3.'+repmat(dy,size(dey3,2),1).*dey3.',2);
                gic1d(indslice,sidx) = nansum(repmat(dx,size(dex1,2),1).*dex1.'+repmat(dy,size(dey1,2),1).*dey1.',2);
            end

            if compute_coupling == 3 %3D
                [avg_strike(indslice,sidx),avg_r(indslice,sidx),std_strike(indslice,sidx),...
                    emag(indslice,sidx),simple_LV(indslice,sidx)] ...
                    = efield_transmissionLine_coupling(exnn3d(s(1:end-1),:),eynn3d(s(1:end-1),:),dx,dy,d);
       
            elseif compute_coupling == 1 %1D 
                [avg_strike(indslice,sidx),avg_r(indslice,sidx),std_strike(indslice,sidx), ...
                    emag(indslice,sidx),simple_LV(indslice,sidx)] ...
                    = efield_transmissionLine_coupling(exnn1d(s(1:end-1),:),eynn1d(s(1:end-1),:),dx,dy,d);
            end


        end
        %disp(['GIC for Transmission Line #',num2str(sidx),' of ',num2str(nlines),' Completed...........................'])
        disp([num2str(100*su_ind(su+1)/length(tind),'%.1f'),'% of Time Steps Completed.............................'])

    end


else %Original option using nested for loops

    disp('LENGTH(TIND)<500. LOOPING OVER EACH TIME STEP-----')

    for sidx = 1:nlines %Loop through lines

        g3d = zeros(length(tind),1);
        g1d = zeros(length(tind),1);

        %Line vertices in latitude (rx) and longitude (ry)
        rx = lines{sidx}(1,:);
        ry = lines{sidx}(2,:);

        centlat = mean(rx);
        centlon = mean(ry);

        %Convert lat long to metres (UTM). Note: y = easting; x = northing 
        % Two methods: Tested both and it doesn't change the results at all.
        [y,x] = geo2utm(ry,rx,ry(round(length(ry)/2)),rx(round(length(rx)/2))); %convert using the transmission line segment midpoint as a reference
        %[y,x] = geo2utm(ry,rx,origlon,origlat); %convert using the full dataset midpoint as an origin (x = northing, y = easting)
        tic

        %Interpolate transmission line values and convert to V/m
        % Here I use scatteredInterpolant to generate 4 functions for ex1d,
        % ey1d, ex3d, and ey3d. Once I have the interpolant function, then I
        % enter the loop and just replace the function values since d.loc never
        % changes. This is faster than e.g. gridding on every time step.

        Fex3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ex3d(tind(1),:).',method);
        Fey3d = scatteredInterpolant(d.loc(:,1),d.loc(:,2),ey3d(tind(1),:).',method);    

        Fex1d = scatteredInterpolant(LAT(:),LON(:),ex1d(tind(1),:).',method);
        Fey1d = scatteredInterpolant(LAT(:),LON(:),ey1d(tind(1),:).',method);

        for tidx = 1:length(tind) %Loop through each time step

            Fex3d.Values = ex3d(tind(tidx),:).';
            exnn3d = Fex3d(rx,ry);

            Fey3d.Values = ey3d(tind(tidx),:).';
            eynn3d = Fey3d(rx,ry);

            Fex1d.Values = ex1d(tind(tidx),:).';
            exnn1d = Fex1d(rx,ry);

            Fey1d.Values = ey1d(tind(tidx),:).';
            eynn1d = Fey1d(rx,ry);


    %         Perform line integral (x and y are in meters, e is in V/m)
    %         Resources:
    %               https://www.mathworks.com/matlabcentral/answers/441416-numerical-calculation-of-line-integral-over-a-vector-field
    %               https://ocw.mit.edu/ans7870/18/18.013a/textbook/HTML/chapter25/section04.html
            g3d(tidx) = (nansum((diff(x)).*(exnn3d(1:end-1)+exnn3d(2:end))/2+(diff(y)).*(eynn3d(1:end-1)+eynn3d(2:end))/2));

            g1d(tidx) = (nansum((diff(x)).*(exnn1d(1:end-1)+exnn1d(2:end))/2+(diff(y)).*(eynn1d(1:end-1)+eynn1d(2:end))/2));

    %         if rem(tidx,3600)==0
    %             disp(['Hour #',num2str(tidx/3600),' Complete'])
    %         end

        end

        gic3d(:,sidx) = g3d;
        gic1d(:,sidx) = g1d;
        toc
        disp(['Transmission Line #',num2str(sidx),' of ',num2str(nlines),' Completed...........................'])

    end

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


        
    