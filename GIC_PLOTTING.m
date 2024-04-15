%% PLOTTING----------------------------------------------------------------
%Plotting routines for Figures 1 through 8 and Supplementary Figures 1
%through 11.

%% Figure 1: Map of zones, transmission lines, cities, and synthetic geoelectric field at 0.01 Hz
%       Also Supplementary Figure 3 and 4 at 0.1 and 0.001 Hz
fidx = nearestpoint(0.01,f);
plot_lim = [47.5 61 -120.5 -109];

plot_zone = true;

%Figure 1a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,101,plot_zone)

if plot_zone
    textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')
end

count = 1; ind = [];
for i = 1:nb
    tmp = find(ininterp(:,i)==1,1);
    if ~isempty(tmp)
        ind(count) = tmp;
        count = count+1;
    end
end
ind = ind(1:14);
ind = indzones;

%refScale = 3000; %good scaling for 0.1 Hz
refScale = 2000; %good scaling for figure at 0.01 Hz
%refScale = 700; %good scaling for figure at 0.001 Hz

quivermc(Clat(indzones)',Clon(indzones)',1000*real(Ey1D(fidx,ind))',1000*real(Ex1D(fidx,ind))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title([num2str(f(fidx)),' Hz'])

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);

cellfun(@(x)plotm(x(1,:),x(2,:),'-k','LineWidth',2), lines);

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
%plotm(49.698,-112.839,'sk','MarkerFaceColor','b','MarkerSize',10); %lethbridge

% Figure 1b: Map of zones, ABT175 and SAB060 locations, and 3-D impedance
initialize_map(plot_lim,zn,provinces,states,102,plot_zone)
%Plot direction of E field for each zone
quivermc(d.loc(:,1),d.loc(:,2),1000*real(Ey3D(fidx,:))',1000*real(Ex3D(fidx,:))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(d.loc(:,1),d.loc(:,2),'.k')

plotm(d.loc(rep,1),d.loc(rep,2),'sk','MarkerSize',10,'MarkerFaceColor','y')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title([num2str(f(fidx)),' Hz'])


%% Figure 2a: Polarization histogram plot
if length(s.rotmag)>1
    figure(201);

    %Calculate angle and radius of each electric field
    th3d = pi/2 - atan2(real(ExR3D),real(EyR3D));
    r3d = sqrt(real(ExR3D).^2+real(EyR3D).^2);

    %Set axis limits based on maximum radius
    upperlim = nanmean(r3d(:))+1*nanstd(r3d(:));
    axlim = ceil(upperlim/(10.^floor(log10(upperlim))))*10.^floor(log10(upperlim));

    %Convert polar angle and radius to cartesian coordinates
    [Y,X] = pol2cart(pi/2-th3d,r3d);

    %Set up 2-D histogram counts
    vXEdge = linspace(-axlim,axlim,50);
    vYEdge = linspace(-axlim,axlim,50);

    mHist2d = hist2d(X(:),Y(:),vYEdge,vXEdge);  

    nXBins = length(vXEdge);
    nYBins = length(vYEdge);
    vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
    vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));


    %Plot 2D Histogram
    pcolor((vXLabel)-0.0001, (vYLabel)-0.0001,log10(mHist2d)); shading flat; hold on
    xlabel('E_y (V/m)');ylabel('E_x (V/m)'); 
    set(gca,'dataaspectratio',[1,1,1])
    c = gray;
    colormap(flip(c(1:round(256/7):end,:))); hcb = colorbar; 
    caxis([-1 2])
    axis([-axlim axlim -axlim axlim])
    grid on
    hcb.Label.String = 'Log_{10} (Count)';


    %Plot average ellipse for the whole dataset
    plot(nanmean(Y,2),nanmean(X,2),'-r','LineWidth',4)

    %Plot ellipses for representative sites ABT175 and SAB060
    plot(nanmean(Y(:,rep(1)),2),nanmean(X(:,rep(1)),2),'-m','LineWidth',2)
    plot(nanmean(Y(:,rep(2)),2),nanmean(X(:,rep(2)),2),'-g','LineWidth',2)
end


%% Figure 2b: Ellipticity map plot of 3-D geoelectric field at 0.01 Hz along with representative sites
%   Representative sites ABT175 (1-D, index 261), SAB060 (3-D, index 476)
%   and ABT272 (Fort Mac, index: 323)
if length(s.rotmag)>1
    ecc = zeros(d.ns,1);
    for i = 1:d.ns
        maxr = max(r3d(:,i));
        minr = min(r3d(:,i));
        ecc(i) = maxr/minr;
    end
    [eccsort, indsort] = sort(ecc);

    plot_lim = [47.5 61 -120.5 -109];
    initialize_map(plot_lim,zn,provinces,states,202,true)

    scatterm(d.loc(indsort,1),d.loc(indsort,2),50*sqrt(eccsort),log10(eccsort),'filled','MarkerEdgeColor','k');
    hcb = colorbar;
    hcb.Label.String = 'Log_{10} (Ellipticity)';
    c = parula;
    colormap(c(1:round(256/7):end,:))
    caxis([0 1.8])
end


%% Figure 3: Map of geoelectric field for real GMD at 0.01 Hz using 1-D and 3-D impedances
%       Also Supplementary Figure 9 and 10 at 0.1 and 0.001 Hz
fidx = nearestpoint(0.010106,f);
%fidx = nearestpoint(0.001,f);
%fidx = 676; %fidx used in paper

refScale = 5; %Scaling for figure at 0.01 Hz in mV/km
%refScale = 0.3; %Scaling for figure at 0.1 Hz
%refScale = 60; %Scaling for figure at 0.001 Hz
refScale = 5;

plot_lim = [47.5 61 -120.5 -109];

%Only plot every 3rd interpolated grid point for visualization purposes
spidx = 1:3:ygrid*xgrid;

%Figure 3a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,301,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

quivermc(LAT(spidx)',LON(spidx)',1000*real(Ey1D(fidx,spidx))',1000*real(Ex1D(fidx,spidx))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(LAT(spidx),LON(spidx),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(['E Field (1-D) @ ',num2str(f(fidx)),' Hz'])

% Figure 3b: 3-D MT Impedance
initialize_map(plot_lim,zn,provinces,states,302,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

quivermc(d.loc(:,1),d.loc(:,2),1000*real(Ey3D(fidx,:))',1000*real(Ex3D(fidx,:))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(d.loc(:,1),d.loc(:,2),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(['E Field (3-D) @ ',num2str(f(fidx)),' Hz'])


%% Figure 4: Histograms of geoelectric field direction in frequency domain for representative MT sites using 1-D and 3-D
fig=figure(400);
clf
set(fig,'Position',[104 463 1536 369])
for i = 1:length(rep)
    is = rep(i);
    is = 167;
    
    %Find nearest interpolated grid point to the reference sites
    [~, indlatlon] = min(distance(d.loc(is,1),d.loc(is,2),LAT(:),LON(:)));

    %Compute the angle between the Ey and Ex component at all frequencies
    th3d = (180/pi)*mod(atan2(real(Ey3D(:,is)),real(Ex3D(:,is))),2*pi);
    th1d = (180/pi)*mod(atan2(real(Ey1D(:,indlatlon)),real(Ex1D(:,indlatlon))),2*pi);

    subplot(1,3,i)
    h(1) = histogram(th1d,360,'FaceColor','g'); hold on
    h(2) = histogram(th3d,360,'FaceColor','m','FaceAlpha',0.3);
    axis([0 360 0 37500])
    legend('1D','3D')
    title(['Frequency domain Angle for Site: ',d.site{is}])
    xlabel('Angle (degrees)')
    ylabel('Count'); grid on
end


%% Figure 5: Map of geoelectric field for real GMD in time domain using 1-D and 3-D impedances
%This time step at 14:02:22 UST has the peak GIC as per the algorithm
%above. In particular, the difference between the 3-D and 1-D models is
%greatest at this time step on the transmission line from Fort McMurray to
%Leismer and also the hypothetical line from Calgary to Nordegg.


tidx = find(b(1).times==datetime('2017-09-08 14:02:22')); %tidx = 28943
tidx = find(b(1).times==datetime('2017-09-08 14:03:12')); %tidx = 28943

%tidx = find(b(1).times==datetime('1989-03-14 01:18:00')); 

%tidx = find(b(1).times==datetime('2021-10-12 10:49:50'));
%tidx = find(b(1).times==datetime('2017-09-08 13:07:03'));

%Plot every 3rd interpolated grid point for visualization purposes
spidx = 1:3:ygrid*xgrid;

%Figure 5a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,501,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

refScale = 0.2; %good scaling for figure at specified time step
refScale = 1000*mean([sqrt(real(ey3d(tidx,:)).^2+ real(ex3d(tidx,:)).^2) sqrt(real(ey1d(tidx,:)).^2+ real(ex1d(tidx,:)).^2)]);
quivermc(LAT(spidx)',LON(spidx)',1000*real(ey1d(tidx,spidx))',1000*real(ex1d(tidx,spidx))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(LAT(spidx),LON(spidx),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(['E Field (1-D) @ ',char(b(1).times(tidx))])

%Figure 5b: 3-D Impedance
initialize_map(plot_lim,zn,provinces,states,502,false)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')
quivermc(d.loc(:,1),d.loc(:,2),1000*real(ey3d(tidx,:))',1000*real(ex3d(tidx,:))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(d.loc(:,1),d.loc(:,2),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(['E Field (3-D) @ ',char(b(1).times(tidx))])

screensize=get(groot,'Screensize');
set(gcf,'Position',[0.5*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])


%% Figure 6: Histograms of geoelectric field direction in time domain for representative MT sites using 1-D and 3-D
fig=figure(600);
clf
%set(fig,'Position',[104 463 1536 369])
for i = 2%:length(rep)
    is = rep(i);
    %is = 183;
    
    %Find nearest interpolated grid point to MT sites
    [~, indlatlon] = min(distance(d.loc(is,1),d.loc(is,2),LAT(:),LON(:)));

    %Compute angle between ey and ex at all time steps
    th3d = (180/pi)*mod(atan2(real(ey3d(:,is)),real(ex3d(:,is))),2*pi);
    th1d = (180/pi)*mod(atan2(real(ey1d(:,indlatlon)),real(ex1d(:,indlatlon))),2*pi);

    %subplot(1,3,i)
    h(1) = histogram(th1d,360,'FaceColor','g'); hold on
    h(2) = histogram(th3d,360,'FaceColor','m','FaceAlpha',0.3);
    axis([0 360 0 13000])
    legend('1D','3D')
    title(['Time Domain domain Angle for Site: ',d.site{is}])
    xlabel('Angle (degrees)')
    ylabel('Count'); grid on
end

%% Figure 7:  Geoelectric time series for representative sites ABT175, SAB060 and ABT272
% This can be used to plot the geoelectric time series for any site by
% changing the "is" variable
is = rep(3); %Use rep(1) for ABT175, rep(2) for SAB060, and rep(3) for ABT272
is = 183;

maxe3d = max(sqrt(ex3d.^2+ey3d.^2),[],1)*1000;

%I plot a slightly more restricted range from 12:00 to 15:00 UT to show the
%peak of the storm more clearly
%tidx = find(b(1).times==datetime('2017-09-08 12:00:00')):find(b(1).times==datetime('2017-09-08 15:00:00')); %12:00 to 15:00 UT  tidx = 21601:32401;
   
%tidx = find(b(1).times==datetime('2017-09-08 13:00:00')):find(b(1).times==datetime('2017-09-08 13:20:00'));
tidx = 1:b(1).nt;

[~,zone_idx] = min(distance(a.loc(is,1),a.loc(is,2),LAT(:),LON(:)));

ccx= corrcoef(ex1d(tidx,zone_idx),ex3d(tidx,is));
ccy= corrcoef(ey1d(tidx,zone_idx),ey3d(tidx,is));

figure(700);
subplot(3,1,1);
plot(b(1).times(tidx),1000*ex1d(tidx,zone_idx),'-b'); hold on
plot(b(1).times(tidx),1000*ex3d(tidx,is),'-r');
ylabel('E (V/km)'); grid on
title(['E_x Time Series for ',d.site{is},'. Correlation Coefficient = ',num2str(ccx(2))]);
legend('1D','3D')

subplot(3,1,2);
plot(b(1).times(tidx),1000*ey1d(tidx,zone_idx),'-b'); hold on
plot(b(1).times(tidx),1000*ey3d(tidx,is),'-r');
ylabel('E (V/km)'); grid on
title(['E_y Time Series for ',d.site{is},'. Correlation Coefficient = ',num2str(ccy(2))])

subplot(3,1,3);
plot(b(1).times(tidx),1000*(ex3d(tidx,is)-ex1d(tidx,zone_idx)),'-k'); hold on
plot(b(1).times(tidx),1000*(ey3d(tidx,is)-ey1d(tidx,zone_idx)),'-g'); hold on
ylabel('E_{3D} - E_{1D} (V/km)')
legend('Diff Ex','Diff Ey');


%% Figure 8: Map of Transmission Line Voltages (Difference between 3-D and 1-D)
plot_lim = [47.5 61 -120.5 -109];

initialize_map(plot_lim,zn,provinces,states,800,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')


tidx = find(b(1).times==datetime('2017-09-08 14:02:00')); 
subindx = tidx-tind(1)+1;
%Set up manual colorbar limits
tmp = dgic(subindx,:);
colormanual(1,:) = [0 0 0.5];
colormanual(2,:) = [0 1 1];
colormanual(3,:) = [0.9 0.9 0];
colormanual(4,:) = [1 0.5 0];
colormanual(5,:) = [1 0 0];
idx = discretize(tmp,[-100 -50 0 50 100 150]);
idx(isnan(idx))=size(colormanual,1);

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

%Loop over lines and plot with the manual color bar
for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat,lon,'-','Color',colormanual(idx(i),:),'LineWidth',2)
end

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
colormap(colormanual)
hcb = colorbar;
hcb.Label.String = 'Difference in Line Voltage (3D - 1D)';
caxis([-100 150]);

%% Figure 9: Forward Modelling Study
% Figures 8a, 8b, and 8c were produced using another set of scripts not
% included here. The models and forward data are included in the forward model
% file. It is relatively easy to load the model and plot slices. 



%% -------------------SUPPLEMENTARY FIGURES--------------------------------
%% Supplementary Figure 1: 1-D models from Trichtchenko et al. (2019)
nbin = length(indzones);
nplots = ceil(sqrt(nbin));
for i = 1:nbin
    subplot(nplots,nplots,i);
    stairs([zn(indzones(i)).rho(1) zn(indzones(i)).rho],[10^-8 zn(indzones(i)).depth(2:end) 10^10]/1000,'k'); hold on
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    axis ij
    grid on
    xticks(10.^(-2:2:6));
    yticks(10.^(-2:1:6));
    axis([0.1 100000 0.1 2000])
    if i>10
        xlabel('Resistivity (\Omega m)')
    end
    if mod(i-1,nplots)==0
        ylabel('Depth (km)')
    end
    title(['Zone ',num2str(indzones(i))])
end

%% Supplementary Figure 2: Impedance data from 1-D models
figure(1); count = 1;
for i = indzones %do 1:11 to reproduce 4.1.3 and 12:nb to reproduce 4.2.3
    c = rand(3,1);
    subplot(2,1,1)
    loglog(d.f,abs(fwd(i).Z),'color',c); hold on
    axis([10^-5 10 10^-4 10^-1]);
    grid on
    ylabel('Impedance Amplitude (Ohms)')
    xlabel('Frequency (Hz)')

    subplot(2,1,2)
    semilogx(d.f,fwd(i).phi,'color',c); hold on
    axis([10^-5 10 0 90]);
    grid on
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    lgstr{count} = ['Zone #',num2str(i)];
    count = count+1;
end

legend(lgstr)

%% Supplementary Figures 3 and 4: See Figure 1 above

%% Supplementary Figure 5: Apparent resistivity and phase curves for ABT175, SAB060 and ABT272
% Figure 5 was produced using another set of scripts not
% included here. The data are still included in the "d" structure and it is
% relatively trivial to plot d.rho and d.pha as a function of d.T.

%% Supplementary Figure 6: Location of Magnetometer Sites

plot_lim = [45 70 -142 -90];
initialize_map(plot_lim,zn,provinces,states,1,false)

plotm([b(1:15).lat],[b(1:15).lon],'.k','MarkerSize',20)
plotm([b([16 17 18 21 22]).lat],[b([16 17 18 21 22]).lon],'.r','MarkerSize',20)
plotm([b([19 20]).lat],[b([19 20]).lon],'.b','MarkerSize',20)
manual_legend('CARISMA','.k','NRCan','.r','USGS','.b');

scale = [zeros(1,19) -2 0 0];
scalelat = [ones(1,20)*0.6 -0.6 0.6];
for i = 1:length(b)
    textm(b(i).lat+scalelat(i),b(i).lon+scale(i),b(i).site,'HorizontalAlignment','center')
end

%% Supplementary Figures 7: Mag time series and spectra for each mag site
%   This function can also be used to view any of the mag sites in the b
%   structure.
tidx = 1:length(b(1).times);

is = 18; %mag site index for Meanook
plot_mag_sites(b,tidx,is,0)

%% Supplementary Figure 8:  Mag Field Frequency Domain
fidx = nearestpoint(0.010106,f);
plot_lim = [47.5 61 -120.5 -109];

spidx = 1:3:ygrid*xgrid;

initialize_map(plot_lim,zn,provinces,states,3,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

refScale = 2000; %good scaling for figure at 0.01 Hz
quivermc(LAT(spidx)',LON(spidx)',real(By_int(spidx,fidx)),real(Bx_int(spidx,fidx)),'color','k','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(LAT(spidx),LON(spidx),'.k')


quivermc(60.5,-118.5,refScale*5,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(refScale*5),' nT'],'HorizontalAlignment','center')
title(['B Field @ ',num2str(f(fidx)),' Hz'])

%% Supplementary Figures 9 and 10: See Figure 3 above

%% Supplementary Table 1:
% excel_table = [gic1d(subindx,:)' gic3d(subindx,:)' dgic(subindx,:)' dgic(subindx,:)'./gic1d(subindx,:)' gtzero'];

%% Supplementary Figure 11: Map of Transmission Line Voltages (1-D)
plot_lim = [47.5 61 -120.5 -109];

%Supplementary Figure 7a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,1,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

tmp = log10(abs([gic1d(subindx,:) gic3d(subindx,:)]));
flr = log10(10);
cel = log10(300);
num = 15;
colormanual = parula(num);
idx = discretize(tmp,flr:(abs(flr-cel)/num):cel);
idx(isnan(idx))=num;

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);


for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat,lon,'-','Color',colormanual(idx(i),:),'LineWidth',3)
    i
end

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
colormap(colormanual)
hcb = colorbar;
hcb.Label.String = 'Log10 (Line Voltage (1D)';
caxis([flr cel])

%% Supplementary Figure 7b: 3-D Impedance
plot_lim = [47.5 61 -120.5 -109];

initialize_map(plot_lim,zn,provinces,states,2,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat,lon,'-','Color',colormanual(idx(i+length(lines)),:),'LineWidth',2)
    i
end

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
%plotm(49.698,-112.839,'sk','MarkerFaceColor','b','MarkerSize',10); %lethbridge
colormap(colormanual)
hcb = colorbar;
hcb.Label.String = 'Log10 (Line Voltage (3D))';
caxis([flr cel])

%% Supplementary Figure 12: Forward Modelling Study (no SABC)
% Supplementary Figure 12 was produced using another set of scripts not
% included here. The models and forward data are included in the forward model
% file. It is relatively easy to load the model and plot slices. 


%% RAW MAGNETOMETER SITES

%Plot raw magnetometer sites
plot_lim = [45 70 -142 -90];
plot_lim = [40 85 -142 -50];
initialize_map(plot_lim,zn,provinces,states,1,false)

plotm([b(:).lat],[b(:).lon],'.k','MarkerSize',20)

scale = [zeros(1,19) -2 0 0];
scalelat = [ones(1,20)*0.6 -0.6 0.6];
for i = 1:length(b)
    textm(b(i).lat+scalelat(i),b(i).lon+scale(i),b(i).site,'HorizontalAlignment','center')
end

%% Plot magnetic field time series at a particular location


bid = 7;
tidx = 1:b(bid).nt;

figure(700);
plot(b(bid).times(tidx),b(bid).x(tidx)-nanmean(b(bid).x(tidx)),'-b'); hold on
plot(b(bid).times(tidx),b(bid).y(tidx)-nanmean(b(bid).y(tidx)),'-r');
ylabel('B (nT)'); grid on
title(['B Time Series for ',b(bid).site])
legend('Bx','By')

%% Plot raw transmission lines
plot_lim = [47.5 61 -120.5 -109];

%Supplementary Figure 7a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,1001,true)
%textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);


for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);

    if i >= 219
        plotm(lat,lon,'-b','LineWidth',3)
    else
        plotm(lat,lon,'-k','LineWidth',2)
    end
    i
    
    %textm(lat(round(length(lat)/2))*1.0002,lon(round(length(lon)/2)),num2str(i),'HorizontalAlignment','center')
end

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river

%% PLOT MAGNETIC FIELD TIME DOMAIN

tidx = find(b(1).times==datetime('2017-09-08 14:02:22')); %tidx = 28943

tidx = find(b(1).times==datetime('1989-03-14 02:08:00')); 

%Plot every 3rd interpolated grid point for visualization purposes
spidx = 1:3:ygrid*xgrid;

%Figure 5a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,501,false)
%textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

refScale = mean(sqrt(bx_int_resh(:).^2+by_int_resh(:).^2)); %good scaling for figure at specified time step
%refScale = 1000*mean([sqrt(real(ey3d(tidx,:)).^2+ real(ex3d(tidx,:)).^2) sqrt(real(ey1d(tidx,:)).^2+ real(ex1d(tidx,:)).^2)]);
quivermc(LAT(spidx)',LON(spidx)',real(by_int_resh(spidx,tidx))',real(bx_int_resh(spidx,tidx))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(LAT(spidx),LON(spidx),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' nT'],'HorizontalAlignment','center')
title(['B Field (1-D) @ ',char(b(1).times(tidx))])

%% Other Figures (Not In Paper): Line Voltage as a function of time for a particular line
% Plot of line voltage as a function of time for a single line

indmax = find(max(abs(gic1d))>100);

linid = indmax(6);

%linid = 219;
linid = 138;

%[~,linid] = ind2sub(size(gic3d),find(abs(gic3d(:))==max(abs(gic3d(:))))); %find line with MAX GIC value
%linid = 228;
%linid=224;

ccx= corrcoef(V1d(:,linid),V3d(:,linid));

tind = 1:b(1).nt;

dl = b(1).times(tind(1));
dr = b(1).times(tind(end));

figure(1000);
%subplot(2,1,1);
plot(b(1).times(tind),V3d(:,linid),'-r'); hold on
plot(b(1).times(tind),V1d(:,linid),'-b');
ylabel('Line Voltage (V)'); grid on
%title(['Line Voltage Time Series for Line ',lineName{linid},'. Correlation Coefficient = ',num2str(ccx(2))],'Interpreter','none');
title(['Line Voltage Time Series for Line ',lineName{linid},'. Max 1-D = ',num2str(max(gic1d(:,linid))),' V/km. Max 3-D = ',num2str(max(gic3d(:,linid))),' V/km'],'Interpreter','none');
datetick('x','HH:MM:ss')
xlabel('Time (MST)')
xlim([dl dr]);
%set(gca,'XLim',[datetime(2012,03,09,0,0,0)-hours(7),datetime(2012,03,09,23,59,59)-hours(7)])
plot([datetime(2012,03,09,2,25,0) datetime(2012,03,09,2,25,0)],[0 max(get(gca,'YLim'))],':k','LineWidth',2)
plot([datetime(2012,03,09,2,25,0) datetime(2012,03,09,2,25,0)],[0 max(get(gca,'YLim'))],':k','LineWidth',2)
manual_legend('1D','-b','3D','-r');
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])

% subplot(2,1,2);
% plot(b(1).times(tind),(gic3d(:,linid)-gic1d(:,linid)),'-k'); hold on
% ylabel('|3D| - |1D| Voltage (V)')
% grid on
% datetick('x','HH:MM:ss')
% xlabel('Time (MST)')
% xlim([dl dr]);
%set(gca,'XLim',[datetime(2012,03,09,0,0,0)-hours(7),datetime(2012,03,09,23,59,59)-hours(7)])

%% Plot maximum line voltage on each line as a function of time

[maxvolt3d, maxvoltind3d] = max(abs(V3d),[],1);
[maxvolt1d, maxvoltind1d] = max(abs(V1d),[],1);

subplot(3,1,3);
plot(b(1).times(maxvoltind3d),maxvolt3d,'o','MarkerFaceColor','b','MarkerEdgeColor','k'); grid on; hold on
plot(b(1).times(maxvoltind1d),maxvolt1d,'v','MarkerFaceColor','r','MarkerEdgeColor','k')
%set(gca,'YScale','log')
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
legend('3D MT','1D Model')
ylabel('Line Voltage (V)')
set(gca,'FontSize',14)
title('Peak Line Voltages')

subplot(3,1,1);
ind = [10 11 18];
for i = ind
    if i == 10
        shiftup = 1500;
    else
        shiftup = 0;
    end
    plot(b(1).times,sqrt(b(i).x.^2+b(i).y.^2)+shiftup,'LineWidth',2); hold on; grid on
end
legend('MCMU','MSTK','MEA')
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
ylabel('Geomagnetic Field (nT)')
title('Geomagnetic Field at 3 Alberta Observatories')
set(gca,'FontSize',14)

subplot(3,1,2);

[~, is(1)] = min(distance(53.5461,-113.4938,d.loc(:,1),d.loc(:,2))); %edmonton
[~, is(2)] = min(distance(51.0477,-114.0719,d.loc(:,1),d.loc(:,2))); %calgary
[~, is(3)] = min(distance(56.726,-111.379,d.loc(:,1),d.loc(:,2))); %fort mac
[~, is(4)] = min(distance(56.234,-117.289,d.loc(:,1),d.loc(:,2))); %peace river

ind = is;
for i = ind
    %plot(b(1).times,1000*mage3d(:,i),'LineWidth',2); hold on;
    plot(b(1).times,1000*nanmean(mage3d,2),'LineWidth',2)
end
grid on
title('Mean Geoelectric Field at all Locations')
%legend('Edmonton','Calgary','Fort McMurray','Peace River')
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
ylabel('Geoelectric Field (V/km)')
%set(gca,'YScale','log')
set(gca,'FontSize',14)

%% Plot peak line voltage as a function of line length
figure; hold on
set(gcf,'Position',[1 58 1088 919])
plot([0 500],[0 500],'--k'); % 1 V/km
plot([0 500],[0 500]*2,'--k'); % 2 V/km
plot([0 500],[0 500]*3,'--k'); % 3 V/km
plot([0 500],[0 500]*4,'--k'); % 4 V/km
plot([0 500],[0 500]*5,'--k'); % 5 V/km
plot([0 500],[0 500]*0.5,'--k'); % 0.5 V/km
plot([0 500],[0 500]*0.1,'--k'); % 0.1 V/km


plot([line_lengths(maxvolt1d>maxvolt3d) line_lengths(maxvolt1d>maxvolt3d)]'/1000,[maxvolt1d(maxvolt1d>maxvolt3d); maxvolt3d(maxvolt1d>maxvolt3d)],'-k'); hold on
plot([line_lengths(maxvolt1d<=maxvolt3d) line_lengths(maxvolt1d<=maxvolt3d)]'/1000,[maxvolt1d(maxvolt1d<=maxvolt3d); maxvolt3d(maxvolt1d<=maxvolt3d)],'-r')
plot(line_lengths/1000,maxvolt3d,'.r','MarkerSize',25); hold on
plot(line_lengths/1000,maxvolt1d,'*k','MarkerSize',10)

ind2plot = [1:length(line_lengths)]'>=219;
ind2plot = line_lengths>250000;
ind2plot = line_lengths>100000 & line_lengths<250000;
ind2plot = line_lengths<100000 & maxvolt3d'>100;
ind2plot = maxvolt3d'>=maxvolt1d';
%ind2plot = maxvolt3d'-maxvolt1d' > 50;
%ind2plot = [];
ind2plot = [1:length(line_lengths)]'==61;

plot([line_lengths(ind2plot) line_lengths(ind2plot)]'/1000,[maxvolt1d(ind2plot); maxvolt3d(ind2plot)],'--m','LineWidth',6); hold on



xlabel('Line Length (km)'); ylabel('Peak Line Voltage (V)')
title(['Lines with Greater Voltage Using MT Impedance = ',num2str(100*sum(maxvolt3d>=maxvolt1d)/length(line_lengths)),'%'])
set(gca,'FontSize',14)
ylim([0 550]);


%Plot corresponding transmission lines based on condition
plot_lim = [47.5 61 -120.5 -109];

%Supplementary Figure 7a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,1001,true)
%textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);

ind = find(ind2plot==1)';
for i = ind
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    
    plotm(lat,lon,'-r','LineWidth',3)
    
    %textm(lat(round(length(lat)/2))*1.0002,lon(round(length(lon)/2)),num2str(i),'HorizontalAlignment','center')
end

plotm(firstpointlat(ind),firstpointlon(ind),'.k','MarkerSize',12)
plotm(lastpointlat(ind),lastpointlon(ind),'.k','MarkerSize',12)
%% Plot Map of Line Voltages larger vs. smaller
plot_lim = [47.5 61 -120.5 -109];

initialize_map(plot_lim,zn,provinces,states,800,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

%Set up manual colorbar limits
tmp = double(maxvolt1d<=maxvolt3d);
colormanual(1,:) = [0 0 0];
colormanual(2,:) = [0 0 0];
colormanual(3,:) = [1 0 0];
colormanual(4,:) = [1 0 0];
colormanual(5,:) = [1 1 1];
idx = discretize(tmp,[-1 0 1 2]);
idx(isnan(idx))=size(colormanual,1);

firstpointlat = cellfun(@(x)x(1,1),lines);
lastpointlat = cellfun(@(x)x(1,end),lines);
firstpointlon = cellfun(@(x)x(2,1),lines);
lastpointlon = cellfun(@(x)x(2,end),lines);



%Loop over lines and plot with the manual color bar
for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat,lon,'-','Color',colormanual(idx(i),:),'LineWidth',3)
end

plotm(firstpointlat,firstpointlon,'.k','MarkerSize',12)
plotm(lastpointlat,lastpointlon,'.k','MarkerSize',12)

manual_legend('Larger Line Voltage using MT impedance','-r','Larger Line Voltage using 1-D Model','-k')
%% Plot outliers: 
% This is all analyzed based on the Sept 8, 2017 storm.
%
%Lines 906 and 928. These lines are the only lines outside of Fort McMurray area
%with both large voltage (>100 V) and greater voltages using MT impedance
%compared to 1-D Models. They also just happen to be in the anomalous area
%with highly polarized geoelectric fields. The geoelectric field magnitude
%is quite small, and yet the 3-D MT impedance gives larger voltages
%   -Note that the 906/928 lines are a double circuit line with termination
% points at Benalto and Sarcee. (Hilariously this is the line that goes
% right past my old house in Calgary along Twelve Mile Coulee Rd...)
%   -The difference in peak voltage is minimal (106 V using MT impedance
% versus 96 V using 1-D model).
%   -% These lines are line indices 143 and 144
%
%Line 12L44: This line is a 500 kV line near Fort Mac and is the only 500
%   kV line with larger voltage using MT impedance. It is index 219.


linid = 219; %line index

subplot(3,1,3)
plot(b(1).times,V1d(:,linid),'-k','LineWidth',2); hold on; grid on
plot(b(1).times,V3d(:,linid),'-r','LineWidth',2);
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
ylabel('Line Voltage (V)')

midpointlat = (firstpointlat(linid)+lastpointlat(linid))/2;
midpointlon = (firstpointlon(linid)+lastpointlon(linid))/2;

[~, is] = min(distance(midpointlat,midpointlon,d.loc(:,1),d.loc(:,2))); % nearest site to the midpoint of Line index iL
[~, ip] = min(distance(midpointlat,midpointlon,LAT(:),LON(:))); % nearest interp grid point to the midpoint of Line index iL

subplot(3,1,1)
plot(b(1).times,1000*mage1d(:,is),'-k','LineWidth',2); hold on;
plot(b(1).times,1000*mage3d(:,ip),'-r','LineWidth',2); 
grid on
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
ylabel('Geoelectric Field (V/km)')
%set(gca,'YScale','log')
set(gca,'FontSize',14)

subplot(3,1,2)
plot(b(1).times,the1d(:,is),'.k','LineWidth',2); hold on;
plot(b(1).times,the3d(:,ip),'.r','LineWidth',2); 
grid on
xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
ylabel('Field Direction (deg)')
%set(gca,'YScale','log')
set(gca,'FontSize',14)
%
%% An idea: Determine the "coupling" between the E-field and the transmission line
%   along it's routing. The idea is to take each line segment and compute
%   the angle between the line segment and the E-field at a given time step
%   If the angle is small (+/-180°), then this suggests strong coupling at
%   that point. Compute the "average coupling" along the line by averaging
%   the angles
%

linid = 219; %12L44 (Thickwood Hills-Livock)
linid = 143; %928L (Benalto-Sarcee)
%linid = 51; %9L15 (Brintnell-Wesley Creek)

tind = 1:length(b(1).times);
compute_coupling = 1;

for linid = 59%1:length(lineName)
    [~,~,avg_strike,avg_r,std_strike,emag,simple_LV] ...
        = calc_line_integral({lines{linid}},tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,'natural',compute_coupling);
    
    [~,indsort] = sort(emag);

    figure(2);
    plot(V3d(:,linid),'-k'); hold on;
    plot(simple_LV,'-r');
    title(['Max P2P LV = ',num2str(max(abs(simple_LV))),' V. Max Integral LV = ',num2str(max(abs(V3d(:,linid)))),' V. % Diff = ',num2str((max(abs(V3d(:,linid)))-max(abs(simple_LV)))/max(abs(V3d(:,linid))))])

    
    close all
    set_figure_size(1);
    if compute_coupling == 1
        scatter(avg_strike(indsort),gic1d(indsort,linid),20,sort(emag)*1000,'filled'); hcb = colorbar;
    elseif compute_coupling == 3
        scatter(avg_strike(indsort),gic3d(indsort,linid),20,sort(emag)*1000,'filled'); hcb = colorbar;
    end
        
    colormap(parula(10)); grid on
    xlabel('Average Angle Between E field and Transmission Line')
    ylabel('Line Voltage (V)')
    hcb.Label.String = 'Geoelectric Field Strength (V/km)';
    set(gca,'FontSize',14)
    %caxis([0 1])
    
    
    [maxe,indmaxe] = max(emag);
    if compute_coupling == 1
        [maxLV(linid),indmaxLV] = max(gic1d(:,linid));
    elseif compute_coupling == 3
        [maxLV(linid),indmaxLV] = max(gic3d(:,linid));
    end
    maxeangle(linid) = avg_strike(indmaxe);
    maxLVangle(linid) = avg_strike(indmaxLV);
    
    potentialMaxLV(linid) = maxe*line_lengths(linid);
    
    
    title(strjoin([lineName{linid}, "Length = ", num2str(line_lengths(linid)/1000), "km.", ...
        "Peak E minus Peak LV) = ", string(b(1).times(indmaxe)-b(1).times(indmaxLV)), ...
        ".Max Peak = ",num2str(potentialMaxLV(linid))," V"]),'interpreter','none');
    
    linestr = lineName{linid};
    linestr = linestr(~ismember(linestr, '/'));
    %print_figure('line_voltage_coupling_1d',linestr)


    plot_lim = [47.5 61 -120.5 -109];

    plot_lim = [min(lines{linid}(1,:))-1 max(lines{linid}(1,:))+1 min(lines{linid}(2,:))-2 max(lines{linid}(2,:))+2]; 

    %Supplementary Figure 7a: 1-D Impedance
    initialize_map(plot_lim,zn,provinces,states,1001,true)
    for i = linid
        lat = lines{i}(1,:);
        lon = lines{i}(2,:);
        plotm(lat,lon,'-','Color','b','LineWidth',3)
    end

    [ex3dint,ey3dint] = interpolate_e_from_Z(d,ex3d,ey3d,LAT,LON);

    tidx = indmaxe;
    %tidx = 50492;
    %Plot every 3rd interpolated grid point for visualization purposes
    spidx = 1:3:ygrid*xgrid;

    inab = find(any(ininterp==1,2));

    spidx = intersect(spidx,inab)';

    for i = 1:length(LAT(:))
        [mintmp(i),indtmp(i)] = min(distance(LAT(i),LON(i),d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84')));

    end


    
    %Figure 5a: 1-D Impedance
    refScale = 1000*mean([sqrt(real(ey3d(tidx,:)).^2+ real(ex3d(tidx,:)).^2) sqrt(real(ey1d(tidx,:)).^2+ real(ex1d(tidx,:)).^2)]);
    quivermc(LAT(spidx)',LON(spidx)',1000*real(ey1d(tidx,spidx))',1000*real(ex1d(tidx,spidx))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
    
    quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
    textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
    
    inclose = find(mintmp<50000);

    spidx = intersect(spidx,inclose);

    %Figure 5b: 3-D Impedance
    quivermc(LAT(spidx)',LON(spidx)',1000*real(ey3dint(tidx,spidx))',1000*real(ex3dint(tidx,spidx))','color','k','reference',refScale,'arrowstyle','tail','linewidth',2);
    
   
    screensize=get(groot,'Screensize');
    set(gcf,'Position',[0.5*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])   

    title(['E Field @ ',char(b(1).times(tidx))])

    manual_legend('1-D Impedance','-r','MT Impedance (interp)','-k')
    
end

%% Another idea
%Calculate the line voltage on every line for 1 V/km E-field oriented at 5
%deg increments. Find the maximum sum total of line voltage to give the
%worst-case direction for the E-field given the network topology

compute_coupling = 0;
tind = 1;
linid = 1:length(lines);

theta = 0:10:360;
V1d = []; V3d = [];
for i = 1:length(theta)

    xprime = 1*cosd(theta(i));
    yprime = 1*sind(theta(i));
    
    ex3d = xprime*ones(1,d.ns);
    ey3d = yprime*ones(1,d.ns);
    
    ex1d = xprime*ones(1,length(LAT(:)));
    ey1d = yprime*ones(size(ex1d));
    
    [V1d(i,:),V3d(i,:)] = calc_line_integral({lines{linid}},tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,'natural',compute_coupling);
end

%% Plot average E-field along transmission line as a function of time

linid = ind(end);
linid = 221;

[~,~,~,~,~,emag1] = calc_line_integral({lines{linid}},tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,'natural',1);
[~,~,~,~,~,emag3] = calc_line_integral({lines{linid}},tind,d,ex3d,ey3d,ex1d,ey1d,LAT,LON,'natural',3);


ccx= corrcoef(emag1,emag3);

tind = 1:b(1).nt;

dl = b(1).times(tind(1));
dr = b(1).times(tind(end));

figure(1000);
plot(b(1).times(tind),emag1*1000,'-b','LineWidth',2); hold on
plot(b(1).times(tind),emag3*1000,'-r');
ylabel('Electric Field Magnitude (V/km)'); grid on
title(['E-field Time Series for Line ',lineName{linid},'. Max 1-D = ',num2str(max(emag1*1000)),' V/km. Max 3-D = ',num2str(max(emag3*1000)),' V/km'],'Interpreter','none');
datetick('x','HH:MM:ss')
xlabel('Time (MST)')
xlim([dl dr]);
%set(gca,'XLim',[datetime(2012,03,09,0,0,0)-hours(7),datetime(2012,03,09,23,59,59)-hours(7)])
%plot([datetime(2012,03,09,2,25,0) datetime(2012,03,09,2,25,0)],[0 max(get(gca,'YLim'))],':k','LineWidth',2)
%plot([datetime(2012,03,09,2,25,0) datetime(2012,03,09,2,25,0)],[0 max(get(gca,'YLim'))],':k','LineWidth',2)
manual_legend('1D','-b','3D','-r');
%xlim([datetime('2017-09-08 12:00:00'),datetime('2017-09-08 15:00:00')])
%% Plot comparison of 1-D vs 3-D Voltages on cross-plot

loglog(maxvolt1d,maxvolt3d,'.k'); hold on; loglog([0.1 1000],[0.1 1000],'--k');
grid on
xlabel('1-D Model Line Voltages (V)')
ylabel('3-D MT Z Line Voltages (V)')
set(gca,'FontSize',14)

figure(2);
plot(maxvolt1d,maxvolt3d,'.k'); hold on; plot([0.1 1000],[0.1 1000],'--k');
grid on
xlabel('1-D Model Line Voltages (V)')
ylabel('3-D MT Z Line Voltages (V)')
set(gca,'FontSize',14)

%% 

%% Plot GIC in Amps along particular line

indmax = find(max(abs(gic3d))>50);

linid = indmax(13);
linid = 220;

dl = b(1).times(tind(1));
dr = b(1).times(tind(end));

figure(1007);
plot(b(1).times(tind),TlineGIC(linid,tind),'-b'); hold on
%plot(b(1).times(tind),V3d(:,linid),'-r');
ylabel('GIC (Amps)'); grid on
title(['GIC Time Series for Line ',lineName{linid}],'Interpreter','none');
datetick('x','HH:MM:ss')
xlabel('Time (MST)')
xlim([dl dr]);

