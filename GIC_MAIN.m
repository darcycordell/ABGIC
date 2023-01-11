%% Set up and User Settings
clc; clearvars;
disp('***********************BEGIN GIC RUN********************************')
full_time = tic;

%---------PRIMARY USER INPUTS TO REPRODUCE FIGURES-------------------------
%To produce Figure 1 (40 - 100 seconds):
%       s.flag = true; %Use synthetic data = true; Use real data = false
%       s.rotmag = 0; %Angle (CW from north) for synthetic B signal (radians)
%       noifft.flag = true; %set as true to *not* perform ifft and only calculate E at a single frequency
%       noifft.freq = 0.001; %frequency in Hz to calculate E
%       
%To produce Figure 2 (~30 minutes):
%       s.flag = true;
%       s.rotmag = linspace(0,2*pi,36); %Angle (CW from north) for synthetic B signal (radians)
%       noifft.flag = true;
%       noifft.freq = 0.01; %frequency in Hz to calculate E
%
%To produce Figure 3 through 8 (17 - 19 minutes):
      s.flag = false;
      s.rotmag = 0;
      noifft.flag = false; 
%
%To Figure 9 is a separate workflow. The models are included in the forward
%modelling folder, but the plotting routines are not included.


%Other user-defined inputs. The defaults are those used in the paper
s.freq = 0.01; % Hz frequency of synthetic signal (0.01 Hz default)
s.mag = 1000; % nT magnitude of synthetic sinusoidal signal (1000 nT default)
s.len = 5000; % Length of synthetic signal (number of samples) (5000 default)

%Mapping and Interpolation Limits
%lim = [45 64 -136 -92]; % lat, lon limits to interpolate
lim = [47.5 61 -120.5 -109];
dlon = 1; % cell sizes in decimal degrees to interpolate
dlat = 1;

%Padding and sample rate for FFT and IFFT used throughout
pad = 10000;
fs = 1; 

if length(s.rotmag)>1 %if you are doing rotations, do not do ifft! (too slow)
    noifft.flag = true;
end

% Load Geographical Info

%Load Shape files with Provincial and State boundaries
%   Located in 01_GEOGRAPHICAL_DATA folder
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

%Load Zones from Trichtchenko et al. (2019)
%   Located in 02_TRICHTCHENKO_ZONES folder
[zn, Clat, Clon, txt, nb] = load_trich_zones;

%>240 kV poswer line network
%   Located in 03_POWERLINE_DATA folder
%lines = load_powerlines('Simplified_Network_36_Lines.kml');

[lines, lineName, line_lengths, voltage, resistance] = post_process_lines;


% Load Time-Domain Mag Data b_s(t) and Convert to Frequency Domain B_s(omega)
%
% Magnetic Field is in nT
%
% "b" is a data structure containing all the information about the magnetic
% fields at the magnetometer site locations. Time-domain fields are lower
% case (e.g. b.x is the time-domain x-component) and frequency-domain
% fields are upper case (e.g. b.X is the frequency-domain x-component).
%
disp('***************LOAD MAGNETIC OBSERVATORY TIME DOMAIN DATA*************')
%Data located in 04_MAG_DATA folder
[b] = load_mag_data(s);

%Delete duplicates
inddelb = find(cellfun(@isempty,({b(:).x})));
b(inddelb) = [];
%
disp('...converting observatory data to frequency domain')
ns = length(b);
%Convert mag sites to the frequency domain.
for i = 1:ns
    [b(i).X,b(i).Y,b(i).f] = calc_fft(b(i).x,b(i).y,fs,pad);
end

f = b(1).f;


% Check for "bad" points
tidx = 1:length(b(1).times);
is = 1;
while 1
    %is = 11; %mag site index for Meanook
    plot_mag_sites(b,tidx,is,0);
    [x,~,click] = ginput(2);
    
    if click~=3
        ax = gca;
        t1 = num2ruler(x(1),ax.XAxis);
        t2 = num2ruler(x(2),ax.XAxis);
        
        [~ ,id1] = min(abs(datenum(t1)-datenum(b(is).times)));
        [~ ,id2] = min(abs(datenum(t2)-datenum(b(is).times)));
        
        b(is).x(id1:id2) = nanmean([b(is).x(1:id1-1);b(is).x(id2+1:end)]);
        b(is).y(id1:id2) = nanmean([b(is).y(1:id1-1);b(is).y(id2+1:end)]);
        b(is).z(id1:id2) = nanmean([b(is).z(1:id1-1);b(is).z(id2+1:end)]);
        
        %b(is).x = inpaint_nans(b(is).x);
        %b(is).y = inpaint_nans(b(is).y);
        %b(is).z = inpaint_nans(b(is).z);
        
        [b(is).X,b(is).Y,b(is).f] = calc_fft(b(is).x,b(is).y,fs,pad);
        
        
    else
        is = is+1;
        
    end
    
    if click(1) == 81 || click(1) == 113
        break
    end
    
    if is>length(b)
        break
    end
    
    close(gcf);
    
end


% Interpolate b_s(t) to b_i(t)
%Time for 1 day recording of 1 Hz data: ~1 minute
disp('*********************INTERPOLATE b(t) SPATIALLY****************')
%Outputs bx_int and by_int(interpolated IAGA B fields)
[bx_int,by_int,LON,LAT] = interpolate_t(b,lim,dlon,dlat);
xgrid = size(bx_int,1);
ygrid = size(bx_int,2);

%
for rotidx = 1:length(s.rotmag)

    if s.flag

        bmag = s.mag*[cos(s.rotmag(rotidx)); sin(s.rotmag(rotidx))]; %set magnitude based on rotation angle
        %Option to include synthetic sinusoidal signal which does not require
        %interpolation
        for i = 1:xgrid
            for j = 1:ygrid
                bx_int(i,j,:) = bmag(1)*cos(2*pi*s.freq*(0:length(b(1).times)-1))';
                by_int(i,j,:) = bmag(2)*cos(2*pi*s.freq*(0:length(b(1).times)-1))';
            end
        end
    end

    %Reshaped for plotting purposes
    by_int_resh = reshape(by_int,size(by_int,1)*size(by_int,2),size(by_int,3));
    bx_int_resh = reshape(bx_int,size(bx_int,1)*size(bx_int,2),size(bx_int,3));



    % CALCULATE B(omega) using Fourier Transform of b_i(t)
    %Time for 672 grid points: ~10 seconds
    disp('******************CONVERT INTERPOLATED b(t) to FREQUENCY-DOMAIN****************')
    nf = b(1).nt+2*pad;
    Bx_int = zeros(xgrid,ygrid,nf); Bx_int = Bx_int+1i*Bx_int;
    By_int = Bx_int;
    for i = 1:xgrid
        for j = 1:ygrid
            [Bx_int(i,j,:),By_int(i,j,:),f] = calc_fft(squeeze(bx_int(i,j,:)),squeeze(by_int(i,j,:)),fs,pad);  
        end
        disp(['Grid Point: (',num2str(i),', ',num2str(j),')'])
    end

    Bx_int = reshape(Bx_int,xgrid*ygrid,nf);
    By_int = reshape(By_int,xgrid*ygrid,nf);


    % Load 3-D MT data from Cedar and 1-D Synthetic Data from Trichtchenko et al. (2019)
    %Takes no time (<0.5 seconds)
    disp('*****************LOAD 1-D AND 3-D MAGNETOTELLURIC IMPEDANCE DATA***************')
    
    %Data located in 05_MT_IMPEDANCE folder
    [d,in,indzones] = load_assign_impedance(zn,'AB_BC_MT_DATA_512_sites.mat');
    d.Z(abs(real(d.Z(:)))>10^5)=NaN;

    [a,zn,fwd] = calc_Z_trich(d,zn,in);

    %Find which interpolated grid points are in each zone polygon
    ininterp = false(length(LAT(:)),nb);
    for i = 1:nb
        ininterp(:,i) = inpolygon(LAT(:),LON(:),zn(i).lat,zn(i).lon);
    end

    %representative sites: ABT175 (1-D), SAB060 (3-D), and ABT272 (Fort Mac)
    rep = [261 475 323]; %Site indices
    %ABT175 = site index 261
    %SAB060 = site index 476
    %ABT272 = site index 323

    disp(['Number of MT sites loaded: ',num2str(d.ns)])
    disp(['Representative sites: ',d.site{rep(1)},', ',d.site{rep(2)},', ',d.site{rep(3)}])

    % CALCULATE E(omega) using 1-D and 3-D methods
    % Frequency domain electric field is in V/m
    disp('********************CALCULATING FREQUENCY-DOMAIN ELECTRIC FIELD****************')
    %This is the biggest calculation of the whole script
    % For the full frequency set and both 1-D and 3-D datasets, this takes
    % around 6 minutes using a parfor loop in calc_E (takes about 20 minutes in serial)
    % Note that I updated calc_E using pagemtimes function and now it only
    % takes 1 minute in serial

    if noifft.flag
        fidx = nearestpoint(noifft.freq,f); %You can do it for only 1 frequency to make
                                     %  it fast, but then you can't do ifft
    else
        fidx = 1:nf;
    end

    tic 
    %Compute E-field using 3-D impedances
    Ex3D = nan(length(fidx),d.ns)+1i*zeros(length(fidx),d.ns); Ey3D = nan(size(Ex3D));
    for i = 1:d.ns
        %Find index of nearest grid point to the site and use that mag field
        %for the calculation
        [~, indlatlon] = min(distance(d.loc(i,1),d.loc(i,2),LAT(:),LON(:)));    
        Zint3D =  interpolate_Z(d.Z(:,:,i),d.f,f);
        [Ex3D(:,i),Ey3D(:,i)] = calc_E(Bx_int(indlatlon,fidx),By_int(indlatlon,fidx),Zint3D(fidx,:));
        disp(['E-field for MT Site #: ',num2str(i),' (Name: ',d.site{i},') completed'])
    end

    if length(s.rotmag)<=1
        %Compute E field using 1-D impedances
        Ex1D = nan(length(fidx),xgrid*ygrid); Ey1D = nan(size(Ex1D));
        for i = 1:xgrid*ygrid

            %Find index of nearest site to the grid point and use that impedance
            %for the calculation
            [~,zone_idx] = min(distance(a.loc(:,1),a.loc(:,2),LAT(i),LON(i)));

            [I,J] = ind2sub(size(LAT),i);

            if ~isempty(zone_idx)      
               Zint1D = interpolate_Z(a.Z(:,:,zone_idx),a.f,f);       
               [Ex1D(:,i),Ey1D(:,i)] = calc_E(Bx_int(i,fidx),By_int(i,fidx),Zint1D(fidx,:));
                disp(['E-field for Grid Point: (',num2str(I),', ',num2str(J),') completed'])       
            end  
        end
    end
    
    if length(s.rotmag)>1
        ExR3D(rotidx,:) = Ex3D;
        EyR3D(rotidx,:) = Ey3D;
    end

    clear Zint3D Zint1D

    toc

end
%%
if ~noifft.flag
    % CALCULATE e(t)
    %
    disp('******************CALCULATING TIME-DOMAIN ELECTRIC FIELD*****************')
    % Time domain electric field is in V/m
    %
    %For the full time series and all 1-D and 3-D locations, this takes
    %approximately 30 seconds
    tic
    ex1d = zeros(b(1).nt,xgrid*ygrid);
    ey1d = zeros(size(ex1d));
    ex3d = zeros(b(1).nt,d.ns);
    ey3d = zeros(size(ex3d));
    for i = 1:xgrid*ygrid
        [ex1d(:,i),ey1d(:,i)] = calc_ifft(Ex1D(:,i),Ey1D(:,i),pad);
        disp(['1-D Impedance E-field for grid point: ',num2str(i)])
    end

    for i = 1:d.ns
        [ex3d(:,i),ey3d(:,i)] = calc_ifft(Ex3D(:,i),Ey3D(:,i),pad);
        disp(['3-D Impedance E-field for MT site: ', num2str(i)])
    end

    mage1d = sqrt(real(ex1d).^2+real(ey1d).^2);
    mage3d = sqrt(real(ex3d).^2+real(ey3d).^2);
    toc

    %% COMPUTE HYPOTHETICAL TRANSMISSION LINE INTEGRAL
    % Get transmission line locations and discretize on a scale smaller than
    % the MT site spacing to capture the geometry of the line if necessary.
    % This will take a VERY long time if you are trying to do it for the entire
    % time series. However, if you just look at a specific window of the time
    % series where the geoelectric field is largest (in our case around 14:00
    % UST), then it is quicker. It takes about 0.3 seconds for each time step
    % so it would take about 7 hours to do all 86400 time steps...but it is
    % only neccesary to do maybe 100 to 1000 time steps to find the peak GIC
    %
    % Sept 28, 2022: Note that I removed griddata function and replaced
    % with scatteredInterpolant where I only interpolate once at the
    % beginning rather than every time step. This speeds things up
    % *considerably* with a single timestep completed in about 0.009
    % seconds. So in this case it would only take about 15 minutes for an
    % entire day of time steps.
    
    %The max geoelectric field occurs throughout the province sometime between
    %13:45 and 14:15 but the peak is spatially distributed e.g. the peak in
    %southern Alberta occurs at a slightly different time than northern
    %Alberta etc. Here, I sweep through the time indices from 28900:29100 and
    %calculate the line integral along the tranmission lines at each time step.

    disp('*******************CALCULATE LINE INTEGRAL VOLTAGE******************')
    lineTic = tic;

    tind = find(b(1).times==datetime('2012-03-09 10:30:00')):find(b(1).times==datetime('2012-03-09 13:30:00')); %13:45 to 14:15 is tind = 27901:29700;
    
    tind = find(b(1).times==datetime('2017-09-08 13:00:00')):find(b(1).times==datetime('2017-09-08 13:20:00')); 
    
    %tind = 1:length(b(1).times);
    
    %tind = find(b(1).times==datetime('2017-09-08 14:02:22')); %tind = 28943; %This is the maximum difference in GIC in restricted range
                                                       %used in the paper
    linid = 1:length(lines);
    %linid = 6:8;
    %linid = 32;
                                                       

    [gic1d,gic3d] = calc_line_integral({lines{linid}},tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,'nearest');
    toc(lineTic);
    %%
    %This results in GIC calculations which have "peak" values at different
    %times depending on the transmission line. Choosing which time step to show
    %is somewhat arbitrary. 
    dgic = (gic3d-gic1d);

    %Find percentage of time steps where 3-D method is larger than 1-D method
    gtzero = sum(dgic>0,1)/length(tind);

    % Peak Difference in GIC magnitude over **restricted interval** (near the
    % peak magnitude)
    [val, ~] = ind2sub(size(dgic(900:end-700,:)),find(dgic(900:end-700,:)==max(max(dgic(900:end-700,:)))));
    subindx = val+900-1; %index for dgic matrix
    tidx = tind(1)+subindx-1; %time index
    b(1).times(tidx) 

    tidx = 28943;
    subindx = tidx+1-tind(1);

    %Supplementary Table 1
    excel_table = [gic1d(subindx,:)' gic3d(subindx,:)' dgic(subindx,:)' dgic(subindx,:)'./gic1d(subindx,:)' gtzero'];

    %%
    for i = 1:length(lines)
        indmax3d(i) = find(gic3d(:,i)==max(gic3d(:,i)));
        indmax1d(i) = find(gic1d(:,i)==max(gic1d(:,i)));
        maxdiffgic(i) = gic3d(indmax3d(i),i)-gic1d(indmax1d(i),i);
    end

end

disp(['TOTAL TIME: ',num2str(toc(full_time))])
disp('*****************************END GIC RUN*****************************')

