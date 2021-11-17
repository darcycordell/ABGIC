%% PLOTTING----------------------------------------------------------------
%Plotting routines for Figures 1 through 8 and Supplementary Figures 1
%through 11.

%% Figure 1: Map of zones, transmission lines, cities, and synthetic geoelectric field at 0.01 Hz
%       Also Supplementary Figure 3 and 4 at 0.1 and 0.001 Hz
fidx = nearestpoint(0.01,f);
plot_lim = [47.5 61 -120.5 -109];

%Figure 1a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,101,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

count = 1; ind = [];
for i = 1:nb
    tmp = find(ininterp(:,i)==1,1);
    if ~isempty(tmp)
        ind(count) = tmp;
        count = count+1;
    end
end

%refScale = 3000; %good scaling for 0.1 Hz
refScale = 2000; %good scaling for figure at 0.01 Hz
%refScale = 700; %good scaling for figure at 0.001 Hz
quivermc(Clat(indzones)',Clon(indzones)',1000*real(Ey1D(fidx,ind))',1000*real(Ex1D(fidx,ind))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title([num2str(f(fidx)),' Hz'])


for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat(1),lon(1),'.k','MarkerSize',12)
    plotm(lat(end),lon(end),'.k','MarkerSize',12)
    plotm(lat,lon,'-k','LineWidth',2)
end

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
%plotm(49.698,-112.839,'sk','MarkerFaceColor','b','MarkerSize',10); %lethbridge

% Figure 1b: Map of zones, ABT175 and SAB060 locations, and 3-D impedance
initialize_map(plot_lim,zn,provinces,states,102,true)
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
    
    %Find nearest interpolated grid point to the reference sites
    [~, indlatlon] = min(distance(d.loc(is,1),d.loc(is,2),LAT(:),LON(:)));

    %Compute the angle between the Ey and Ex component at all frequencies
    th3d = (180/pi)*mod(atan2(real(Ey3D(:,is)),real(Ex3D(:,is))),2*pi);
    th1d = (180/pi)*mod(atan2(real(Ey1D(:,indlatlon)),real(Ex1D(:,indlatlon))),2*pi);

    subplot(1,3,i)
    h(1) = histogram(th1d,36,'FaceColor','g'); hold on
    h(2) = histogram(th3d,36,'FaceColor','m','FaceAlpha',0.3);
    axis([0 360 0 7500])
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
tidx = 28943;

%Plot every 3rd interpolated grid point for visualization purposes
spidx = 1:3:ygrid*xgrid;

%Figure 5a: 1-D Impedance
initialize_map(plot_lim,zn,provinces,states,501,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

refScale = 0.2; %good scaling for figure at specified time step
%refScale = 1000*mean([sqrt(real(ey3d(tidx,:)).^2+ real(ex3d(tidx,:)).^2) sqrt(real(ey1d(tidx,:)).^2+ real(ex1d(tidx,:)).^2)]);
quivermc(LAT(spidx)',LON(spidx)',1000*real(ey1d(tidx,spidx))',1000*real(ex1d(tidx,spidx))','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
plotm(LAT(spidx),LON(spidx),'.k')

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(['E Field (1-D) @ ',char(b(1).times(tidx))])

%Figure 5b: 3-D Impedance
initialize_map(plot_lim,zn,provinces,states,502,true)
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
set(fig,'Position',[104 463 1536 369])
for i = 1:length(rep)
    is = rep(i);
    
    %Find nearest interpolated grid point to MT sites
    [~, indlatlon] = min(distance(d.loc(is,1),d.loc(is,2),LAT(:),LON(:)));

    %Compute angle between ey and ex at all time steps
    th3d = (180/pi)*mod(atan2(real(ey3d(:,is)),real(ex3d(:,is))),2*pi);
    th1d = (180/pi)*mod(atan2(real(ey1d(:,indlatlon)),real(ex1d(:,indlatlon))),2*pi);

    subplot(1,3,i)
    h(1) = histogram(th1d,36,'FaceColor','g'); hold on
    h(2) = histogram(th3d,36,'FaceColor','m','FaceAlpha',0.3);
    axis([0 360 0 12000])
    legend('1D','3D')
    title(['Time Domain domain Angle for Site: ',d.site{is}])
    xlabel('Angle (degrees)')
    ylabel('Count'); grid on
end

%% Figure 7:  Geoelectric time series for representative sites ABT175, SAB060 and ABT272
% This can be used to plot the geoelectric time series for any site by
% changing the "is" variable
is = rep(3); %Use rep(1) for ABT175, rep(2) for SAB060, and rep(3) for ABT272

%I plot a slightly more restricted range from 12:00 to 15:00 UT to show the
%peak of the storm more clearly
tidx = 21601:32401;

ccx= corrcoef(ex1d(tidx,is),ex3d(tidx,is));
ccy= corrcoef(ey1d(tidx,is),ey3d(tidx,is));

figure(700);
subplot(3,1,1);
plot(b(1).times(tidx),1000*ex1d(tidx,is),'-b'); hold on
plot(b(1).times(tidx),1000*ex3d(tidx,is),'-r');
ylabel('E (V/km)'); grid on
title(['E_x Time Series for ',d.site{is},'. Correlation Coefficient = ',num2str(ccx(2))]);
legend('1D','3D')

subplot(3,1,2);
plot(b(1).times(tidx),1000*ey1d(tidx,is),'-b'); hold on
plot(b(1).times(tidx),1000*ey3d(tidx,is),'-r');
ylabel('E (V/km)'); grid on
title(['E_y Time Series for ',d.site{is},'. Correlation Coefficient = ',num2str(ccy(2))])

subplot(3,1,3);
plot(b(1).times(tidx),1000*(ex3d(tidx,is)-ex1d(tidx,is)),'-k'); hold on
ylabel('E_{3D} - E_{1D} (V/km)')


%% Figure 8: Map of Transmission Line Voltages (Difference between 3-D and 1-D)
plot_lim = [47.5 61 -120.5 -109];

initialize_map(plot_lim,zn,provinces,states,800,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

%Set up manual colorbar limits
tmp = dgic(subindx,:);
colormanual(1,:) = [0 0 0.5];
colormanual(2,:) = [0 1 1];
colormanual(3,:) = [0.9 0.9 0];
colormanual(4,:) = [1 0.5 0];
colormanual(5,:) = [1 0 0];
idx = discretize(tmp,[-100 -50 0 50 100 150]);
idx(isnan(idx))=size(colormanual,1);

%Loop over lines and plot with the manual color bar
for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat(1),lon(1),'.k','MarkerSize',12)
    plotm(lat(end),lon(end),'.k','MarkerSize',12)
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

tmp = log10([gic1d(subindx,:) gic3d(subindx,:)]);
flr = log10(10);
cel = log10(300);
num = 15;
colormanual = jet(num);
idx = discretize(tmp,flr:(abs(flr-cel)/num):cel);
idx(isnan(idx))=num;

for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat(1),lon(1),'.k','MarkerSize',12)
    plotm(lat(end),lon(end),'.k','MarkerSize',12)
    plotm(lat,lon,'-','Color',colormanual(idx(i),:),'LineWidth',2)
end

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
plotm(56.234,-117.289,'sk','MarkerFaceColor','b','MarkerSize',10); %peace river
colormap(colormanual)
hcb = colorbar;
hcb.Label.String = 'Log10 (Line Voltage (1D)';
caxis([flr cel])

%Supplementary Figure 7b: 3-D Impedance
plot_lim = [47.5 61 -120.5 -109];

initialize_map(plot_lim,zn,provinces,states,2,true)
textm(Clat*1.005,Clon,txt,'HorizontalAlignment','center')

for i = 1:length(lines)
    lat = lines{i}(1,:);
    lon = lines{i}(2,:);
    plotm(lat(1),lon(1),'.k','MarkerSize',12)
    plotm(lat(end),lon(end),'.k','MarkerSize',12)
    plotm(lat,lon,'-','Color',colormanual(idx(i+length(lines)),:),'LineWidth',2)
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

