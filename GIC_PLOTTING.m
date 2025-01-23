
%% Base Map
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

dbase = load('AB_BC_MT_DATA_526_sites_230321.mat');

plot_lim = [47.5 61 -120.5 -109];

screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',14)

geoshow(states,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5); hold on
geoshow(provinces,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5)

plotm(dbase.d.loc(:,1),dbase.d.loc(:,2),'.','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5)

gs = scaleruler;
setm(gs,'FontSize',16,'MajorTick',[0 100 200])
for i = 1:nLines
    lat = L(i).Loc(:,1);
    lon = L(i).Loc(:,2);
    if L(i).Voltage == 500
        tcolor = 'b';
    else
        tcolor = 'k';
    end
    plotm(lat,lon,'-','Color',tcolor,'LineWidth',2)
end

subLoc = reshape([S(:).Loc],2,length(S))';
plotm(subLoc(:,1),subLoc(:,2),'o','MarkerFaceColor','w','MarkerEdgeColor','k'); 


plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac

%% SUPPLEMENTARY FIGURE 1: MAGNETIC OBSERVATORY LOCATIONS
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

plot_lim = [39 70 -142 -90];

screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',14)

geoshow(states,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5); hold on
geoshow(provinces,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5)


for i = 1:length(b)-1
    textm(b(i).lat,b(i).lon+1*sind(b(i).lat),b(i).site,'FontSize',14,'FontWeight','bold')
end

textm(b(end).lat,b(end).lon-5*sind(b(end).lat),b(end).site,'FontSize',14,'FontWeight','bold')


plotm([b(1:9).lat],[b(1:9).lon],'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10)
plotm([b(10:14).lat],[b(10:14).lon],'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
plotm([b(15:end).lat],[b(15:end).lon],'o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10)
setm(gca,'FontSize',14);


%% FIGURE 1: MAP OF STUDY AREA
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

plot_lim = [47.5 61 -120.5 -109];

screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',14)

geoshow(states,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5); hold on
geoshow(provinces,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5)

plotm(d.loc(:,1),d.loc(:,2),'v','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5)



for i = 1:nLines
    lat = L(i).Loc(:,1);
    lon = L(i).Loc(:,2);

    if L(i).Voltage>400
        plotm(lat,lon,'-','Color','b','LineWidth',3)
    else
        plotm(lat,lon,'-','Color','k','LineWidth',3)
    end
end

plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac

subLoc = reshape([S(:).Loc],2,length(S))';
plotm(subLoc(:,1),subLoc(:,2),'o','MarkerFaceColor','w','MarkerEdgeColor','k'); 

isub(1) = find(strcmp({S.Name},'ELLERSLIE 89S'));
isub(2) = find(strcmp({S.Name},'KEEPHILLS 320P'));
isub(3) = find(strcmp({S.Name},'SUNNYBROOK 510S'));
isub(4) = find(strcmp({S.Name},'CROSSINGS 511S'));
isub(5) = find(strcmp({S.Name},'BENNETT 520S'));

for i = 1:length(isub)
    plotm(subLoc(isub(i),1),subLoc(isub(i),2),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',12); 
end


plotm([b(:).lat],[b(:).lon],'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',10)



setm(gca,'FontSize',14)
%% SUPPLEMENTARY FIGURE 8A: Map of Geoelectric Field at Peak Time of Bowmanton

%exu = xmax*ones(1,d.ns);
%eyu = ymax*ones(1,d.ns);

sitelist = 0; count = 1; siteind = [];
for i = 1:d.ns
    if ~any(sitelist==i)
        siteind(count) = i;
        count = count+1;
        D = distance(d.loc(i,1),d.loc(i,2),d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84'));

        sitelist = [sitelist; find(D<25000 & D~=0)];
    end
end

%siteind = 1:d.ns;


indsub = min(find(contains({S(:).Name},'HORIZON UP')));
[~,ind] = max(abs(GIC_Subs(indsub,:)));

%ind = indmax;

ind = find(b(1).times==GIC_Times(ind));
%[~,ind] = min(abs(datetime(2003,10,29,6,48,0)-b(1).times));
%ind = 95443;
%ind = 341;
%ind = find(b(1).times==GIC_Times(itime));

%ind = find(b(1).times==datetime(2023,4,24,7,19,50));

%ind = find(b(1).times==datetime(2023,4,24,2,15,39));
ind = find(b(1).times==datetime(2024,5,11,2,18,0));

%ind = find(b(1).times == datetime(2012,3,9,9,20,00));
refScale = 1/10;

for i = 1%maxeGroups(15) %1:1:length(ind)
    openfig('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Network_Model_Scripts/basemap.fig');
    setm(gca,'FontSize',16)
    quivermc(d.loc(siteind,1),d.loc(siteind,2),1000*ey(ind(i),siteind)',1000*ex(ind(i),siteind)','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
    %quivermc(d.loc(siteind,1),d.loc(siteind,2),1000*ymax(ind(i),siteind)',1000*xmax(ind(i),siteind)','color','k','reference',refScale,'arrowstyle','tail','linewidth',2);
    %quivermc(d.loc(siteind,1),d.loc(siteind,2),1000*ey1(ind(i),siteind)',1000*ex1(ind(i),siteind)','color','k','reference',refScale,'arrowstyle','tail','linewidth',2);
    

    quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
    textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
    title(char(b(1).times(ind(i))-tz,'yyyy-MMM-dd HH:mm:ss'));

%    plotm(xlat,xlon,'oy')
   
    
end
set(gca,'FontSize',16)
%cd(runPath);
%% FIGURE 8B: Geoelectric Field Direction and Magnitude Plot

screensize=get(groot,'Screensize');
is = 269;
eMag = sqrt(ex(:,is).^2+ey(:,is).^2);
eTheta = (180/pi)*mod(atan2(ey(:,is),ex(:,is)),2*pi);

datelimits = [t1 t2];

fig = figure(251);
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.5*screensize(3) 0.3*screensize(4)])
[~,indsort] = sort(eMag,'ascend');
scatter(b(1).times(indsort),eTheta(indsort),500000*eMag(indsort),1000*eMag(indsort),'filled'); hold on
%caxis([min(datenum(b(1).times)) max(datenum(b(1).times))])
colormap(parula(5))
hcb = colorbar; hcb.Label.String = 'Electric Field Magnitude (V/km)';
ylim([0 360])
xlim(datelimits)
ylabel('Electric Field Bearing (°E of N)')
grid on; box on
set(gca,'FontSize',16)

[maxebow,indmaxbow] = max(eMag);

%title(['Max E field of ',num2str(maxebow),' V/km @ ',datestr(b(1).times(indmaxbow)),' with bearing ',num2str(eTheta(indmaxbow)),'° E of N'])


%% FIGURE 2: Simplified DC-equivalent network circuit 
% AESO PDF SLD


%% FIGURE 3-7 GEOMAGNETIC FIELD, GEOELECTRIC FIELD, GIC MODEL AND HALL 
halltrans = 6;

tightOn = 0;
plotdelta = 0;


%t1 = datetime('2021-10-12 02:00:00');
%t2 = datetime('2021-10-12 16:00:00');

%t1 = datetime('2023-04-24 01:00:00');
%t2 = datetime('2023-04-24 09:00:00');


%t1 = datetime('2023-03-23 06:00:00');
%t2 = datetime('2023-03-23 18:00:00');

%t1 = datetime('2024-05-11 8:00:00');
%t2 = datetime('2024-05-11 11:00:00');


%t1 = GIC_Times(1);
%t2 = GIC_Times(end);
%t2 = datetime(2024,9,17,6,0,0);

tz = hours(0);

datelimits = [t1 t2];
nHall = size(GIC_Hall,2);
tshift = 0;
%GIC_Times = GIC_Times_orig + seconds(tshift);
tgic = tgic_orig - seconds(tshift);
%No corrections except up-down bias
if contains(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T1') || ...
    contains(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T2')
    hall = (-(GIC_Hall{:,2}+GIC_Hall{:,3}))-mean(-(GIC_Hall{:,2}+GIC_Hall{:,3}),'omitnan');
    %hall = -GIC_Hall{:,2}-mean(-GIC_Hall{:,2});
    %hall2 = -GIC_Hall{:,3}-mean(-GIC_Hall{:,3});
    tran1 = -GIC_Hall{:,2};
    tran1([49758 49771 49780 49785 49800 49805 49808 49814 49830 49834 49842 49799 49804 49807 49813 49818 49803 49798 49819 49824 49786 49781 49820 49821 49835 49843]) = NaN;
    tran1 = inpaint_nans(tran1);

    tran2 = -GIC_Hall{:,3};
    tran2([49768 49778 49796 49804 49809 49819 49820 49822 49769 49777 49797 49805 49808 49807 49810]) = NaN;
    tran2 = inpaint_nans(tran2);
    hall = tran1+tran2-mean(tran1+tran2,'omitnan');
    

elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x12S_T1')
    hall = -GIC_Hall{:,halltrans}-mean(-GIC_Hall{:,halltrans},'omitnan');
    
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x320P_T6')
    hall = -1*GIC_Hall{:,halltrans}-mean(-GIC_Hall{:,halltrans},'omitnan');

elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x510S_T11')
    hall = -1*GIC_Hall{:,halltrans}-mean(-GIC_Hall{:,halltrans},'omitnan');

elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x511S_T1')
    hall = -1*GIC_Hall{:,halltrans}-mean(-GIC_Hall{:,halltrans},'omitnan');

elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x520S_T1')
    hall = -1*GIC_Hall{:,halltrans}-mean(-GIC_Hall{:,halltrans},'omitnan');

end

%High pass GIC data doesn't chance much.
%hall = highpass(hall,1/10000,1/2);

%1000 Ωm halfspace
% if strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T1') || ...
%     strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T2')
%     hall = (-(GIC_Hall{:,2}+GIC_Hall{:,3})+4);
%     GIC_Times = GIC_Times_orig + seconds(60);
%     %hall = (-(GIC_Hall{:,2})+4);
%     %hall = -(GIC_Hall{:,3})+4;
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x12S_T1')
%     hall = -GIC_Hall{:,halltrans};
%     GIC_Times = GIC_Times_orig + seconds(0);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x320P_T6')
%     hall = -6*GIC_Hall{:,halltrans}-4;
%     GIC_Times = GIC_Times_orig + seconds(60);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x510S_T11')
%     hall = -3*GIC_Hall{:,halltrans}+12;
%     GIC_Times = GIC_Times_orig + seconds(0);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x511S_T1')
%     hall = -2*GIC_Hall{:,halltrans}+6;
%     GIC_Times = GIC_Times_orig + seconds(60);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x520S_T1')
%     hall = -3*GIC_Hall{:,halltrans}-7;
%     GIC_Times = GIC_Times_orig + seconds(60);
% end

%MT impedance
% if strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T1') || ...
%     strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T2')
%     hall = -0.5*(GIC_Hall{:,2}+GIC_Hall{:,3})+4;
%     GIC_Times = GIC_Times_orig + seconds(60);
%     %hall = (-(GIC_Hall{:,2})+4);
%     %hall = -(GIC_Hall{:,3})+4;
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x12S_T1')
%     hall = -GIC_Hall{:,halltrans};
%     GIC_Times = GIC_Times_orig + seconds(0);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x320P_T6')
%     hall = -2*GIC_Hall{:,halltrans}-0;
%     GIC_Times = GIC_Times_orig + seconds(60);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x510S_T11')
%     hall = -1*GIC_Hall{:,halltrans}+5;
%     GIC_Times = GIC_Times_orig + seconds(0);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x511S_T1')
%     hall = -1*GIC_Hall{:,halltrans}+5;
%     GIC_Times = GIC_Times_orig + seconds(0);
% elseif strcmp(GIC_Hall.Properties.VariableNames{halltrans},'x520S_T1')
%     hall = -1*GIC_Hall{:,halltrans}-0;
%     GIC_Times = GIC_Times_orig + seconds(60);
% end

%plot(tgic,-(GIC_Hall{:,3}),'-r','LineWidth',2); hold on
%plot(tgic,-(GIC_Hall{:,2}),'-b','LineWidth',2); hold on

if contains(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T1') %ellerslie 89S
    isub = find(strcmp({S.Name},'ELLERSLIE 89S'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x89S_T2') %ellerslie 89S
    isub = find(strcmp({S.Name},'ELLERSLIE 89S'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x12S_T1') %heartland 12S
    isub = find(strcmp({S.Name},'HEARTLAND 12S'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x320P_T6') %keephills 320P
    isub = find(strcmp({S.Name},'KEEPHILLS 320P'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x510S_T11') %sunnybrook 510S
    isub = find(strcmp({S.Name},'SUNNYBROOK 510S'));
    %isub = find(strcmp({S.Name},'GENESEE 330P'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x511S_T1') %crossings 511S
    isub = find(strcmp({S.Name},'CROSSINGS 511S'));
elseif contains(GIC_Hall.Properties.VariableNames{halltrans},'x520S_T1') %bennett 520S
    isub = find(strcmp({S.Name},'BENNETT 520S'));
end

%isub = find(contains({S.Name},'HORIZON UP')); isub = isub(1);
% [maxvalnew,itimemaxnew] = max(abs(Snew(isub).GIC));
% 3*maxvalnew
[maxval,itimemax] = max(abs(S(isub).GIC));
% 3*maxval
% 
plot_substation(GIC_Times-tz,S,L,T,GIC_Subs,GIC_Lines,GIC_Trans,isub,itimemax,1)
% plot_substation(GIC_Times-tz,Snew,Lnew,Tnew,g,l,t,isub,itimemax,1)


[~,ia,ib] = intersect(GIC_Times,tgic);

%Uncomment to not plot GIC data and only plot model
%ia = []; ib = [];

if isempty(ia)
%if 1
    tgic = GIC_Times;
    hall = nan(size(S(isub).GIC))';
    [~,ia,ib] = intersect(GIC_Times,tgic);
end


tbw1 = find(isbetween(GIC_Times(ia),t1,t2));
tbw2 = find(isbetween(tgic(ib),t1,t2));

%hall = hall/3;

g = squeeze(3*S(isub).GIC(ia(tbw1)));
h = hall(ib(tbw2));
%h2 = -GIC_Hall{ib(tbw2),5}-mean(-GIC_Hall{:,5});
%g = squeeze(3*s1d(isub,ia(tbw1)));
%g = squeeze(3*su(isub,ia(tbw1)));
%g = squeeze(3*GIC_Trans(84,1,ia(tbw1)));%+squeeze(3*GIC_Trans(154,1,ia(tbw1)));
%g2 = squeeze(3*S(86).GIC(ia(tbw1)));
%g = gicG(ia(tbw1))';

% figure(100);
% plot3(h,g,1:length(g),'.k'); hold on; grid on; view([0 90])
% plot3(h,g1d,1:length(g),'.r'); hold on; grid on; view([0 90])
% 
% plot(xlim,xlim,'--r')
% cc = corrcoef(GIC_Subs(isub,ia(tbw1)),s1d(isub,ia(tbw1)));
% set(gca,'FontSize',16);
% xlabel('GIC (A/phase) [Measured]')
% ylabel('GIC (A/phase) [Modelled]')
% legend('3-D','1-D')
% title([S(isub).Name,' (R = ',num2str(cc(2)),', \beta = ',num2str(rc(isub)),')'])

%hall_model_space(:,halltrans) = g;
%ex_model_space(:,halltrans) = exSub(:,isub);
%ey_model_space(:,halltrans) = eySub(:,isub);
%g = squeeze(3*GIC_Trans(isub,2,ia(tbw1)));%+squeeze(3*GIC_Trans(154,2,ia(tbw1)));

stats = statistics(hall(ib(tbw2)),g',ones(length(tbw1),1));
%stats = statistics(hall(ib(tbw2)),g1d',ones(length(tbw1),1));
cc= corrcoef(hall(ib(tbw2)),g);
%cc = corrcoef(hall(ib(tbw2)),g1d);
P = 1 - stats.rms/std(hall(ib(tbw2)),'omitnan');

% Compute the mean of your data
yhat = mean(g);
%yhat = mean(g1d);
xhat = mean(hall(ib(tbw2)));
% Compute regression coefficients in the least-square sense
beta = (hall(ib(tbw2))' - xhat)*(g' - yhat)/sum((hall(ib(tbw2)) - xhat).^2); % Regression coefficient
%beta = (hall(ib(tbw2))' - xhat)*(g1d' - yhat)/sum((hall(ib(tbw2)) - xhat).^2); % Regression coefficient


[maxGICdata,indmaxData] = max(abs(hall(ib(tbw2))));
[maxGICmodel,indmaxModel] = max(abs(g));
%[maxGICmodel,indmaxModel] = max(abs(g1d));
maxGICdataTime = tgic(ib(tbw2(indmaxData)));
maxGICmodelTime = GIC_Times(ia(tbw1(indmaxModel)));
perdiff = (maxGICmodel' - maxGICdata)/maxGICdata;

%
% hallC = []; counti = 1; countj = 1;
% for i = [2 5:nHall]
%     for j = [2 5:nHall]
%         cchall = corrcoef(GIC_Hall{:,i},GIC_Hall{:,j});
% 
%         
%         hallC(counti,countj) = cchall(2);
% 
%         countj = countj+1;
%     end
%     counti = counti+1;
%     countj = 1;
% end
% 
% hallC = triu(hallC,1);
% 
% ind = find(hallC~=0);
% 
% [I,J] = ind2sub(size(hallC),ind);
% 
% %plot(J,hallC(ind),'*k'); grid on
% 
%     hall1 = (-(GIC_Hall{:,2}+GIC_Hall{:,3}))-mean(-(GIC_Hall{:,2}+GIC_Hall{:,3}));
%     %hall = -GIC_Hall{:,2}-mean(-GIC_Hall{:,2});
%     %GIC_Times = GIC_Times_orig + seconds(60);
% 
%     hall2 = -1*GIC_Hall{:,5}-mean(-GIC_Hall{:,5});

%plot(tgic,-hall1/2); hold on
%plot(tgic,hall2);
%xlim(datelimits)
%

fontSize = 18;

figure;
set(gcf,'Position',[680 59 987 918]);
if plotdelta
    ax1 = subplot(4,1,3);
else
    ax1 = subplot(3,1,3);
end
plot(tgic(ib(tbw2))-tz,hall(ib(tbw2)),'-r','LineWidth',2); hold on
plot(GIC_Times(ia(tbw1))-tz,g,'-k','LineWidth',0.5); hold on %neutral

%plot(dmm.time,dmm.I_y,'-r','LineWidth',2);  grid on; hold on;
%plot(GIC_Times,gdmm,'-k','LineWidth',2);

%plot(GIC_Times(ia(tbw1)),3*su(isub,ia(tbw1)),'--b','LineWidth',1.5); hold on %neutral
%plot(GIC_Times(ia(tbw1)),3*s1d(isub,ia(tbw1)),':g','LineWidth',1.5); hold on %neutral
%plot(GIC_Times(ia(tbw1)),g2,'-b','LineWidth',1.5); hold on %neutral
%plot(GIC_Times(ia(tbw1)),hall_model_space(:,halltrans),'-m')
%plot(GIC_Times,3*squeeze(GIC_Trans(isub,1,:)),'-b','LineWidth',2); hold on %HV (Series) Winding
%plot(GIC_Times,3*squeeze(GIC_Trans(isub,2,:)),'--g','LineWidth',2); hold on %LV(Common) Winding

%plot(tgic(ib(tbw2)),hall2(ib(tbw2)),'-b','LineWidth',2);


yl = ylim;

%plot([GIC_Times(ia(tbw1(indmaxModel))) GIC_Times(ia(tbw1(indmaxModel)))],yl,'--k')
%plot([tgic(ib(tbw2(indmaxData))) tgic(ib(tbw2(indmaxData)))],yl,'--r')
grid on
ylabel('GIC (A)')
legend('Hall Probe','Neutral Modelled')
%legend('Neutral Modelled');
%legend('Hall Probe','Data-Space Method','Model-Space Method')
%legend('Hall Probe','3-D Model','1-D Piecewise','AutoUpdate','off')
%legend('T1','T2')
%legend('Sub Neutral','HV (Series)', 'LV (Common)','Hall Probe')
%text(ax1,0.05,0.9,'C','FontSize',30,'FontWeight','Bold','Units','normalized')
set(gca,'FontSize',fontSize)
xlim(datelimits-tz);
% ax = gca;
% ax.XTick = t1:minutes(1):t2;
% ax.XTickLabelMode = 'manual';
% ax.XTickLabel = datestr(ax.XTick,'HH:MM'); % see datestr formatOut
%set(gca,'Xticklabel',[])
if tightOn
    set(gca,'Xticklabel',[])
    
else

    %title(['Total GIC For Substation: ',S(isub).Name,'. RMS = ',num2str(stats.rms),'. r = ',num2str(cc(2)),'. P = ',num2str(P),'. \beta = ',num2str(beta)])
    title(['Total GIC For Substation: ',S(isub).Name,'. RMS = ',num2str(stats.rms),'. R = ',num2str(cc(2)),'. \beta = ',num2str(beta)])
    
end


% subplot(2,1,2);
% plot(GIC_Times(ia(tbw1)),hall(ib(tbw2))-g','-k','LineWidth',2); hold on %difference
% grid on
% ylabel('\Delta GIC Measured - Predicted (A)')
% %legend('Sub Neutral','HV (Series)', 'LV (Common)','Hall Probe')
% set(gca,'FontSize',fontSize)
% xlim(datelimits);
if plotdelta
    ax2 = subplot(4,1,1);
else
    ax2 = subplot(3,1,1);
end

[~,indquerye] = min(distance(subLoc(isub,1),subLoc(isub,2),d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84')));
[~,indquery] = min(distance(subLoc(isub,1),subLoc(isub,2),[b(:).lat],[b(:).lon],referenceEllipsoid('wgs84')));
indquery = 3;


plot(b(1).times-tz,b(indquery).x-mean(b(indquery).x,'omitnan'),'-r','LineWidth',1.5); hold on
plot(b(1).times-tz,b(indquery).y-mean(b(indquery).y,'omitnan'),'-b','LineWidth',1.5); grid on
%plot(b(1).times(1:end-1),diff(bxSub(isub,:)),'-r','LineWidth',1.5); hold on
%plot(b(1).times(1:end-1),diff(bySub(isub,:)),'-b','LineWidth',1.5); grid on
ylabel('B (nT)')
legend('B_x','B_y','AutoUpdate','off')
set(gca,'FontSize',fontSize)
%text(ax2,0.05,0.9,'A','FontSize',30,'FontWeight','Bold','Units','normalized')
xlim(datelimits-tz);
ylim([-2600 1300])
%title(S(isub).Name)
%text(b(indquery).site)
% ax = gca;
% ax.XTick = t1:minutes(1):t2;
% ax.XTickLabelMode = 'manual';
% ax.XTickLabel = datestr(ax.XTick,'HH:MM'); % see datestr formatOut

set(gca,'Xticklabel',[])

if tightOn
    
else

    title(['Max GIC Data = ',num2str(maxGICdata),' A. Max GIC Model = ',num2str(maxGICmodel),' A. Percent Diff = ',num2str(100*perdiff),'%'])
end

if plotdelta
    ax3 = subplot(4,1,2);
else
    ax3 = subplot(3,1,2);
end
%plot(b(1).times-tz,1000*ex(:,indquerye),'-r','LineWidth',1.5); hold on
%plot(b(1).times-tz,1000*ey(:,indquerye),'-b','LineWidth',1.5); grid on; grid(gca,'minor')

plot(b(1).times-tz,1000*mage(:,indquerye),'-k','LineWidth',1.5); grid on; grid(gca,'minor'); hold on


%plot(b(1).times,1000*mage(:,indquerye-10),'-r','LineWidth',1.5);
%plot(b(1).times,1000*ex_model_space(:,halltrans),'-g','LineWidth',1.5);
%plot(b(1).times,1000*ey_model_space(:,halltrans),'-m','LineWidth',1.5);

ylabel('|E| (V/km)')
%legend('E_x','E_y','AutoUpdate','off')

%legend('E_x (Data Space)','E_y (Data Space)','E_x (Model Space)','E_y (Model Space)')
set(gca,'FontSize',fontSize)
%text(ax3,0.05,0.9,'B','FontSize',30,'FontWeight','Bold','Units','normalized')
xlim(datelimits-tz);
ylim([0 1.3])
% ax = gca;
% ax.XTick = t1:minutes(1):t2;
% ax.XTickLabelMode = 'manual';
% ax.XTickLabel = datestr(ax.XTick,'HH:MM'); % see datestr formatOut
set(gca,'Xticklabel',[])
if tightOn
    %set(gca,'Xticklabel',[])
else
    %title('Conductivity Model: MT impedance')
end

if plotdelta
    ax4 = subplot(4,1,4);
    plot(tgic(ib(tbw2))-tz,hall(ib(tbw2))-g','-b','LineWidth',1.5); hold on
    yl = ylim;
    % ax = gca;
    % ax.XTick = t1:minutes(1):t2;
    % ax.XTickLabelMode = 'manual';
    % ax.XTickLabel = datestr(ax.XTick,'HH:MM'); % see datestr formatOut
    
    grid on
    ylabel('\Delta GIC (A)')
    set(gca,'FontSize',fontSize)
    text(ax4,0.05,0.9,'D','FontSize',30,'FontWeight','Bold','Units','normalized')
    xlim(datelimits-tz);
else
    ax4 = ax1;
end

if tightOn
    ax2.Position = [ax2.Position(1:3) 0.25];
    ax1.Position = [ax1.Position(1:3) 0.25];
    ax3.Position = [ax3.Position(1:3) 0.25];
    ax4.Position = [ax4.Position(1:3) 0.25];
end

%ax4.XTickLabelMode = 'manual';
%ax4.XTickLabel = datestr(ax4.XTick,'HH:MM');

% if tz == hours(0)
%     text(ax4,0.8,-0.25,['UT ',datestr(ax4.XTick(end),'mmm dd, yyyy')],'FontSize',18,'Units','normalized')
% else
%     text(ax4,0.8,-0.25,['MLT ',datestr(ax4.XTick(end),'mmm dd, yyyy')],'FontSize',18,'Units','normalized')
% end



% subplot(3,1,1);
% ylim1 = ylim;
% plot([datetime(2024,5,11,8+6,16,15) datetime(2024,5,11,8+6,16,15)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,10+6,8,52) datetime(2024,5,11,10+6,8,52)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,15+6,2,42) datetime(2024,5,11,15+6,2,42)],ylim,'--k','LineWidth',2);
% ylim(ylim1);
% 
% subplot(3,1,2);
% ylim1 = ylim;
% plot([datetime(2024,5,11,8+6,16,15) datetime(2024,5,11,8+6,16,15)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,10+6,8,52) datetime(2024,5,11,10+6,8,52)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,15+6,2,42) datetime(2024,5,11,15+6,2,42)],ylim,'--k','LineWidth',2);
% ylim(ylim1);
% 
% subplot(3,1,3);
% ylim1 = ylim;
% plot([datetime(2024,5,11,8+6,16,15) datetime(2024,5,11,8+6,16,15)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,10+6,8,52) datetime(2024,5,11,10+6,8,52)],ylim,'--k','LineWidth',2);
% plot([datetime(2024,5,11,15+6,2,42) datetime(2024,5,11,15+6,2,42)],ylim,'--k','LineWidth',2);
% ylim(ylim1);

linkaxes([ax1 ax2 ax3 ax4],'x')
%% SAVE AUTOMATICALLY FIGURES
% folder = '/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Network_Model_Scripts/Results/20230423/Revision';
% filename = ['SFig02_',strrep(S(isub).Name,' ','_'),'_MT_ModEM_noDiff'];
% printfile=[folder,'/',filename];
% print('-depsc','-painters',[printfile,'.eps'])

%% SUPPLEMENTARY FIGURE 8: CROSS PLOTS

figure(97)
subplot(3,2,3)
plot(hall(ib(tbw2)),g,'.k'); hold on
plot(xlim,ylim,'-r')
axis equal; grid on
xlabel('GIC Measured (A)');
ylabel('GIC Modelled (A)')
title(S(isub).Name)
set(gca,'FontSize',14)

%% SUPPLEMENTARY FIGURE 8: HIST PLOT (DID NOT USE)
figure(98)
subplot(3,2,4)
[N,Xedge,Yedge] = histcounts2(hall(ib(tbw2)),g',-200:200,-200:200);
pcolor(Xedge(1:end-1),Yedge(1:end-1),log10(N)'); hold on
%plot(hall(ib(tbw2)),g','.k')
plot(xlim,ylim,'-r')
axis equal; grid on; shading flat
xlabel('GIC Measured (A)');
ylabel('GIC Modelled (A)')
title(S(isub).Name)
xlim([min(hall(ib(tbw2))) max(hall(ib(tbw2)))])
ylim([min(g) max(g)])
set(gca,'layer','top')
set(gca,'FontSize',14)
colorbar
cmap = flip(gray(5));
caxis([-0.5 2]);
colormap(cmap(1:end,:))


%% FIGURE 8A: MAP OF PEAK GIC AT EACH SUBSTATION 
[maxSub,indmaxSub] = max(abs(GIC_Subs),[],2);
for i = 1:nSubs
    maxnon(i) = 3*GIC_Subs(i,indmaxSub(i));
end
maxnon(maxnon==0) = NaN;
maxSub(maxSub==0) = NaN;

plot_lim = [47.5 61 -120.5 -109];
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);
plot_zone = true;


%Figure 1a: 1-D Impedance
screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',18)

geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)
geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)

for i = 1:nLines
    lat = L(i).Loc(:,1);
    lon = L(i).Loc(:,2);
    plotm(lat,lon,'-','Color','k','LineWidth',3); hold on
end

%indpeak = 35641;
%indnotzero = abs(TneutGIC(:,indpeak))>0;

subLoc = reshape([S(:).Loc],2,length(S))';

[sortmax,indsort] = sort(maxSub,'ascend');
othermax = 3*maxSub;
othermax(isnan(othermax))=0;
[sortmax,indsort2] = sort(othermax,'descend');

maxSub = 3*maxSub; %did not mulitply by three for 2024 paper whoops
scatterScale = 10; %50 for paper

scatterm(subLoc(indsort,1),subLoc(indsort,2),maxSub(indsort)*scatterScale,(maxnon(indsort))','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
%textm(subLoc(indsort2(1:10),1),subLoc(indsort2(1:10),2),num2cell(1:10))

bmk = [150 100 50];
bmk = [50 25 5];
scatterm([59 58 57],[-119 -119 -119],bmk*scatterScale,bmk','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
textm([59 58 57],[-119 -119 -119],num2cell(bmk),'FontSize',16,'HorizontalAlignment','center')
%h = scatterm(subLoc(:,1),subLoc(:,2),3*abs(GIC_Subs(:,indmaxModel)*10)+10,3*(GIC_Subs(:,indmaxModel)'),'filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%scatterm(subLoc(indnotzero,1),subLoc(indnotzero,2),abs(TneutGIC(indnotzero,indpeak))*4,TneutGIC(indnotzero,indpeak)','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%title(['Peak Neutral GIC'])
colorlims = [-100 100];
caxis(colorlims)
colormap(redblue(diff(colorlims)*2/10,colorlims))
%caxis([-0.1 0.1])
set(gca,'FontSize',18)

%% FIGURE 8B: Plots of peak GIC and mag observatories through time
fig = figure(17);
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) 0.65*screensize(4)])

scatter(GIC_Times(indmaxSub)-tz,othermax,150,subLoc(:,1),'filled','MarkerEdgeColor','k'); hold on
%text(GIC_Times(indmaxSub(indsort2([1:10])))-tz,othermax(indsort2([1:10])),num2cell([1:10]));
text(GIC_Times(indmaxSub(indsort2([1:10])))-tz+hours(0.5),othermax(indsort2([1:10])),{S(indsort2([1:10])).Name});
caxis([50 60]); 
colormap(gray(5))
hcb =colorbar('Location','southoutside'); hcb.Label.String = 'Substation Latitude';
xlim(datelimits-tz); grid on
set(gca,'FontSize',18)
set(gca,'layer','bottom')
ylabel('Total Peak GIC (A)')
box on

fig = figure(18);
set(fig,'Position',[0.1*screensize(3) 1*screensize(4) 0.35*screensize(3) 0.25*screensize(4)])

plot(b(1).times-tz,sqrt(b(5).x.^2+b(5).y.^2)-mean(sqrt(b(5).x.^2+b(5).y.^2)),'-b','LineWidth',3); hold on
plot(b(1).times-tz,sqrt(b(6).x.^2+b(6).y.^2)-mean(sqrt(b(6).x.^2+b(6).y.^2)),'-k','LineWidth',3);hold on
plot(b(1).times-tz,sqrt(b(9).x.^2+b(9).y.^2)-mean(sqrt(b(9).x.^2+b(9).y.^2)),'-r','LineWidth',3);
legend('MCMU','MSTK','VULC')
xlim(datelimits-tz); grid on
set(gca,'XTickLabels',[])
set(gca,'FontSize',18)
ylabel('Magnetic Field (nT)')


%%
magsitelist = {'MCMU','MSTK','VULC'};
for i = 1:3
    ax(i) = subplot(3,1,i);
    imag = find(contains({b(:).site},magsitelist{i}));
    plot(b(1).times-tz,sqrt(b(imag).x.^2+b(imag).y.^2)-mean(sqrt(b(imag).x.^2+b(imag).y.^2)),'-k','LineWidth',3); hold on
    %plot(b(1).times-tz,b(imag).x-mean(sqrt(b(imag).x)),'-r','LineWidth',3); hold on
    %plot(b(1).times-tz,b(imag).y-mean(sqrt(b(imag).y)),'-b','LineWidth',3); hold on

    title(b(imag).site)
    xlim(datelimits-tz); grid on
    set(gca,'FontSize',18)

    if i == 2
    ylabel('Horizontal Magnetic Field (nT)')
    end

end

linkaxes(ax)


%% FIGURE XX: PLOTS OF GEOELECTRIC FIELD AND GIC AT PEAK TIMES


f1 = openfig('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Network_Model_Scripts/basemap_noMT_subplot.fig');
set(f1,'Position',[271 539 714 432]);
set(f1,'Position',[1 50 1264 930]);

subplot(1,2,1)
refScale = 0.1; %good scaling for 1 V/km E field
quivermc(60.3,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.7,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center','FontSize',12)

subplot(1,2,2)
bmk = [50 25 5];
scatterm([59 58 57],[-119 -119 -119],bmk*50/3,bmk','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
textm([59 58 57],[-119 -119 -119]+2,num2cell(bmk),'FontSize',16,'HorizontalAlignment','center')

indsub = find(contains({S(:).Name},'HORIZON UPGRAD'),1);
%indsub = find(contains({S(:).Name},'BOWMANTON'));
%indsub = find(contains({S(:).Name},'HORIZON'),1);
%indsub = find(contains({S(:).Name},'KEEPHILLS'),1);
%indsub = find(contains({S(:).Name},'THICKWOOD'),1);


[~,indg] = max(abs(GIC_Subs(indsub,:)));
ind = find(b(1).times==GIC_Times(indg));
    
    [~,indsort] = sort(abs(GIC_Subs(:,indg)),'ascend');

    vecplot = 3*GIC_Subs(indsort,indg);
    vecplot(vecplot==0) = NaN;
    
    subplot(1,2,2)
    h1 = scatterm(subLoc(indsort,1),subLoc(indsort,2),abs(vecplot*50/3),(vecplot)','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Substation GIC (A)';
    colorlims = [-70 70];
    caxis(colorlims)
    colormap(redblue(diff(colorlims)*2/10,colorlims))
    ax = gca;
    ax.Position = [0.5703 0.1100 0.3347 0.8150];
    
    setm(gca,'FontSize',18)
    set(gca,'FontSize',18)
    drawnow

    subplot(1,2,1)
    h2 = quivermc(d.loc(:,1),d.loc(:,2),1000*ey(ind,:)',1000*ex(ind,:)','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
    %title(char(b(1).times(ind),'yyyy-MMM-dd HH:mm:ss'));
    setm(gca,'FontSize',18)
    set(gca,'FontSize',18)

    sgtitle([char(b(1).times(ind)-tz,'yyyy-MMM-dd HH:mm:ss')],'FontSize',18','FontWeight','Bold');
    drawnow
    


%% SUPPLEMENTARY VIDEO

%Base map subplots
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

dbase = load('AB_BC_MT_DATA_526_sites_230321.mat');

plot_lim = [47.5 61 -120.5 -109];

screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.7*screensize(3) screensize(4)])

for iplot = 1:2
    subplot(1,2,iplot)
    worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
    setm(gca,'FontSize',14)
    
    geoshow(states,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5); hold on
    geoshow(provinces,'DefaultFaceColor',[0.9 1 0.7],'DefaultEdgeColor','black','facealpha',0.5)
    
    %plotm(dbase.d.loc(:,1),dbase.d.loc(:,2),'.','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5)
    
    for i = 1:nLines
        lat = L(i).Loc(:,1);
        lon = L(i).Loc(:,2);
        if L(i).Voltage == 500
            tcolor = 'b';
        else
            tcolor = 'k';
        end
        plotm(lat,lon,'-','Color',tcolor,'LineWidth',2)
    end
    
    subLoc = reshape([S(:).Loc],2,length(S))';
    plotm(subLoc(:,1),subLoc(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','k'); 
    
    
    plotm(53.5461,-113.4938,'sk','MarkerFaceColor','b','MarkerSize',10) %edmonton
    plotm(51.0477,-114.0719,'sk','MarkerFaceColor','b','MarkerSize',10) %calgary
    plotm(56.726, -111.379,'sk','MarkerFaceColor','b','MarkerSize',10); %fort mac
end
    
indsub = find(contains({S(:).Name},'KEEPHILLS'),1,'first');
[~,indg] = max(abs(GIC_Subs(indsub,:)));
%indg = 1;
%indg = find(GIC_Times=='2023-04-24 07:26:46');
%indg = find(GIC_Times=='2024-05-11 09:40:00');
%indg = find(GIC_Times=='2012-03-09 09:20:00');
%indg = find(GIC_Times=='2023-03-24 03:15:46');
%indg = find(GIC_Times=='2022-02-03 11:47:45');
ind = find(b(1).times==GIC_Times(indg));
subLoc = reshape([S(:).Loc],2,length(S))';

vidFile = VideoWriter('test.mp4','MPEG-4');
vidFile.FrameRate = 16;
open(vidFile);

f1 = openfig('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Network_Model_Scripts/basemap_noMT_subplot.fig');
%scatterm([59 60 61],[-118.5 -118.5 -118.5],[10 30 50]*50,([10 30 50])','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Substation GIC (A)';

subplot(1,2,1)
refScale = 0.1; %good scaling for 1 V/km E field
quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')

subplot(1,2,2)
bmk = [100 50 25];
bmk = [60 30 15];
scatterScale = 30;
scatterm([59 58 57],[-119 -119 -119],bmk*scatterScale,bmk','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
textm([59 58 57],[-119 -119 -119]+0.75,num2cell(bmk),'FontSize',16,'HorizontalAlignment','center')

iloop = [max([1 indg-60*10]):1:min([size(GIC_Subs,2) indg+60*10])] - indg;
k = 1;
%for i = -2000:2:2000
%for i = -300:300
for i = iloop

    
    [~,indsort] = sort(abs(GIC_Subs(:,indg+i)),'ascend');

    vecplot = 3*GIC_Subs(indsort,indg+i);
    vecplot(vecplot==0) = NaN;
    
    subplot(1,2,2)
    h1 = scatterm(subLoc(indsort,1),subLoc(indsort,2),abs(vecplot*scatterScale),(vecplot)','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Substation GIC (A)';
    colorlims = [-70 70];
    caxis(colorlims)
    colormap(redblue(diff(colorlims)*2/10,colorlims))
    ax = gca;
    ax.Position = [0.5703 0.1100 0.3347 0.8150];
    title(char(b(1).times(i+ind)-tz,'yyyy-MMM-dd HH:mm:ss'));
    drawnow

    subplot(1,2,1)
    itime = i+ind;
    h2 = quivermc(d.loc(:,1),d.loc(:,2),1000*ey(itime,:)',1000*ex(itime,:)','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
    title(char(b(1).times(itime)-tz,'yyyy-MMM-dd HH:mm:ss'));
    drawnow
    
    F(k) = getframe(gcf);
    delete(h1); delete(h2);
    writeVideo(vidFile,F(k));
    k = k+1;
end
close(vidFile);





%% PLOT OF PEAK GEOELECTRIC FIELDS IN MAP

%% FIGURE 8A: MAP OF PEAK GIC AT EACH SUBSTATION 
[maxSub,indmaxSub] = max(mage,[],1);
for i = 1:d.ns
    maxnon(i) = 1000*mage(indmaxSub(i),i);
end
maxnon(maxnon==0) = NaN;
maxSub(maxSub==0) = NaN;

plot_lim = [47.5 61 -120.5 -109];
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);
plot_zone = true;


%Figure 1a: 1-D Impedance
screensize=get(groot,'Screensize');
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',18)

geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)
geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)

for i = 1:nLines
    lat = L(i).Loc(:,1);
    lon = L(i).Loc(:,2);
    plotm(lat,lon,'-','Color','k','LineWidth',3); hold on
end

%indpeak = 35641;
%indnotzero = abs(TneutGIC(:,indpeak))>0;


[sortmax,indsort] = sort(maxSub,'ascend');
othermax = 3*maxSub;
othermax(isnan(othermax))=0;
[sortmax,indsort2] = sort(othermax,'descend');

maxSub = 1000*maxSub; %did not mulitply by three for 2024 paper whoops
scatterScale = 100; %50 for paper

scatterm(d.loc(indsort,1),d.loc(indsort,2),maxSub(indsort)*scatterScale,(maxnon(indsort))','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
%textm(subLoc(indsort2(1:10),1),subLoc(indsort2(1:10),2),num2cell(1:10))

bmk = [150 100 50];
bmk = [15 30 60];
scatterm([59 58 57],[-119 -119 -119],bmk*scatterScale,bmk','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
textm([59 58 57],[-119 -119 -119],num2cell(bmk),'FontSize',16,'HorizontalAlignment','center')
%h = scatterm(subLoc(:,1),subLoc(:,2),3*abs(GIC_Subs(:,indmaxModel)*10)+10,3*(GIC_Subs(:,indmaxModel)'),'filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%scatterm(subLoc(indnotzero,1),subLoc(indnotzero,2),abs(TneutGIC(indnotzero,indpeak))*4,TneutGIC(indnotzero,indpeak)','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%title(['Peak Neutral GIC'])
colorlims = [0 6];
caxis(colorlims)
colormap(redblue(7,colorlims))
%caxis([-0.1 0.1])
set(gca,'FontSize',18)





%%
f = figure(1);
screensize=get(groot,'Screensize');
set(f,'Position',[0.1*screensize(3) 0*screensize(4) 0.8*screensize(3) screensize(4)])

for i = 1%:d.ns
    subplot(2,1,1)
    plot(b(1).times,1000*ex(:,i),'-r','LineWidth',1.5); hold on
    plot(b(1).times,1000*ey(:,i),'-b','LineWidth',1.5);
    ylabel('E Field (V/km)')
    grid on
    set(gca,'FontSize',16);
    legend('E_x','E_y');
    title(['Site #',sprintf('%03d',i),': ',d.site{i}])
    
    subplot(2,1,2)
    plot(b(1).times,bx(i,:)-mean(bx(i,:),'omitnan'),'-r','LineWidth',1.5); hold on
    plot(b(1).times,by(i,:)-mean(by(i,:),'omitnan'),'-b','LineWidth',1.5);
    ylabel('B Field (nT)')
    grid on
    set(gca,'FontSize',16);
    legend('B_x','B_y');

%     if plot_fig
%         exportgraphics(f,['site',sprintf('%03d',i),'_',d.site{i},'.png'],'Resolution',100)
%         clf(f);
%     end
    
    

end


%%














plot(GIC_Subs(isub,:),'-k'); hold on; 

%Series winding transformers are identical current
plot(squeeze(GIC_Trans(isub,1,:)),'-b'); 
plot(squeeze(GIC_Trans(154,1,:)),'--r');

%Common winding transformers are identical current
plot(squeeze(GIC_Trans(isub,2,:)),'-g');
plot(squeeze(GIC_Trans(154,2,:)),'--m');


%%
%GIC neutral to ground at substation is identical to the sum of the common
%winding transformers
plot(GIC_Subs(isub,:),'-k'); hold on; 
plot(squeeze(GIC_Trans(isub,2,:))+squeeze(GIC_Trans(154,2,:)),'--g');


%%
plot(squeeze(GIC_Subs(isub,:))-squeeze(GIC_Trans(isub,2,:))','-k'); hold on; 

%% Geoelectric Field Direction and GIC Magnitude Plot

screensize=get(groot,'Screensize');
[~,is] = min(distance(subLoc(isub,1),subLoc(isub,2),d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84')));
eMag = sqrt(ex(:,is).^2+ey(:,is).^2);
eTheta = (180/pi)*mod(atan2(ey(:,is),ex(:,is)),2*pi);

[~,ia,ib] = intersect(b(1).times,tgic);

datelimits = [t1 t2];

fig = figure(251);
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.5*screensize(3) 0.3*screensize(4)])
[~,indsort] = sort(hall(ib),'ascend');
scatter(b(1).times(ia),eTheta(ia),5*abs(hall(ib)),hall(ib),'filled'); hold on
%caxis([min(datenum(b(1).times)) max(datenum(b(1).times))])
colormap(parula(5))
hcb = colorbar; hcb.Label.String = 'Electric Field Magnitude (V/km)';
ylim([0 360])
xlim(datelimits)
ylabel('Electric Field Bearing (°E of N)')
grid on; box on
set(gca,'FontSize',16)

[maxebow,indmaxbow] = max(eMag);


%% Plot polar diagrams E theta vs GIC etc
[vald,indsort] = sort(distance(subLoc(isub,1),subLoc(isub,2),d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84')),'ascend');

isub1 = find(strcmp({S.Name},'ELLERSLIE 89S'));
isub2 = find(strcmp({S.Name},'KEEPHILLS 320P'));

is = indsort(1);
eMag = sqrt(ex(:,is).^2+ey(:,is).^2);
eTheta = (180/pi)*mod(atan2(ey(:,is),ex(:,is)),2*pi);

[~,ia,ib] = intersect(b(1).times,tgic);

val = abs(hall(ib));
val2 = abs(hall2(ib));
%val(val<0) = 0;

val3 = 1000*meanE(ia);
%val2(val2<0) = 0;
figure;
set(gcf,'Position',[1 48 984 929])
subplot(2,2,3)
%polarplot(eTheta(ia)*pi/180,val,'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
polarplot(avg_strike(ia)*pi/180,val,'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
%polarplot(avg_strike(ia)*pi/180,abs(S(isub).GIC(ib)),'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
%polarplot(eTheta(ia)*pi/180,abs(S(isub).GIC(ib)),'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';

title([S(isub1).Name,' GIC (A)'])
set(gca,'FontSize',14)

subplot(2,2,4)
%polarplot(eTheta(ia)*pi/180,val,'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
polarplot(avg_strike(ia)*pi/180,val2,'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
%polarplot(avg_strike(ia)*pi/180,abs(S(isub).GIC(ib)),'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
%polarplot(eTheta(ia)*pi/180,abs(S(isub).GIC(ib)),'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';

title([S(isub2).Name,' GIC (A)'])
set(gca,'FontSize',14)


subplot(2,2,1)
%polarplot(eTheta(ia)*pi/180,1000*eMag(ia),'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
polarplot(avg_strike(ia)*pi/180,val3,'.k'); pax = gca; pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
title('Mean Alberta Geoelectric Field (V/km)')
set(gca,'FontSize',14)

subplot(2,2,2)
rose_geog(-avg_strike(ia),30,5000,'r')
%rose_geog(-eTheta(ia),30,5000,'b')
title('Geoelectric Field Azimuth Histogram')
set(gca,'FontSize',14)

%title(['Max E field of ',num2str(maxebow),' V/km @ ',datestr(b(1).times(indmaxbow)),' with bearing ',num2str(eTheta(indmaxbow)),'° E of N'])


%% Plot DMM

curdir = pwd;
cd('/Users/darcycordell/Library/CloudStorage/SynologyDrive-Mac/Work/Projects/GIC/000_DATA/DMM/20230424')
dmm = load('20230424_DMM.mat');
dmm.time.TimeZone = '';
cd(curdir);
%Compare GIC in Line to DMM measurement

dmm.I_y = highpass(dmm.I_y,1/10000,1);

[~,ia,ib] = intersect(GIC_Times,dmm.time);
%[~,ia,ib] = intersect(tgic,dmm.time);

tbw1 = find(isbetween(GIC_Times(ia),t1,t2));
%tbw1 = find(isbetween(tgic(ia),t1,t2));
tbw2 = find(isbetween(dmm.time(ib),t1,t2));

iline1 = find(contains({L.Name},'1206'));
iline2 = find(contains({L.Name},'1212'));
g = 3*L(iline1).GIC+3*L(iline2).GIC;
%g = 3*lu(225,:)+3*lu(225,:);
%g = hall;
g_stats = g(ia(tbw1));

maxGICdata = max(abs(dmm.I_y(ib(tbw2))));
maxGICmodel = max(abs(g));
perdiff = (maxGICmodel - maxGICdata)/maxGICdata;

stats = statistics(dmm.I_y(ib(tbw2)),g_stats',ones(length(tbw1),1));
cc= corrcoef(dmm.I_y(ib(tbw2)),g_stats,'rows','complete');
P = 1 - stats.rms/std(dmm.I_y(ib(tbw2)),'omitnan');

% Compute the mean of your data
yhat = mean(g_stats);
xhat = mean(dmm.I_y(ib(tbw2)));
% Compute regression coefficients in the least-square sense
beta = (dmm.I_y(ib(tbw2))' - xhat)*(g_stats' - yhat)/sum((dmm.I_y(ib(tbw2)) - xhat).^2); % Regression coefficient


figure(1);
set(gcf,'Position',[680 59 987 918]);
%ax4 = gca;
ax4 = subplot(3,1,3);
datelimits = [t1 t2];
%plot(tgic,g,'-k','LineWidth',2); hold on;
plot(dmm.time,dmm.I_y,'-r','LineWidth',2); grid on; hold on
plot(GIC_Times,g,'-k','LineWidth',0.5); hold on;

xlim(datelimits)
ylabel('GIC (A)')
%ylim([-80 120]);


legend('DMM','Line GIC Modelled')
%title(['Total GIC For Lines: ',L(225).Name,' and ',L(226).Name,'. RMS = ',num2str(stats.rms),'. r = ',num2str(cc(2)),'. P = ',num2str(P),'. \beta = ',num2str(beta)])
title(['Total GIC For Lines: ',L(225).Name,' and ',L(226).Name,'. RMS = ',num2str(stats.rms),'. R = ',num2str(cc(2)),'. \beta = ',num2str(beta)])
%text(ax4,0.05,0.9,'C','FontSize',30,'FontWeight','Bold','Units','normalized')
%text(t1+hours(1),-45,['Corr Coef = ',num2str(cc(2)),'. R.M.S. = ',num2str(stats.rms)],'FontSize',18)
set(gca,'FontSize',18)
ax4.XTickLabel = datestr(ax4.XTick,'HH:MM');
text(ax4,0.8,-0.25,['UT ',datestr(ax4.XTick(end),'mmm dd, yyyy')],'FontSize',18,'Units','normalized')


%subplot(3,1,1);
%title(['Max GIC Data = ',num2str(maxGICdata),' A. Max GIC Model = ',num2str(maxGICmodel),' A. Percent Diff = ',num2str(100*perdiff),'%'])


%% FIGURE XX FOR 2022 02 03 STORM PAPER
%Plot geoelectric field and GIC at 5 time steps

t1 = datetime('2022-02-03 11:43:31');
t2 = datetime('2022-02-03 11:47:41');
t3 = datetime('2022-02-03 11:52:11');
t4 = datetime('2022-02-03 11:55:16');
t5 = datetime('2022-02-03 11:57:36');

ind = find(b(1).times==t5);

refScale = 0.1; %good scaling for 1 V/km E field

openfig('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Network_Model_Scripts/basemap_noMT_subplot.fig');

subplot(1,2,1)
quivermc(d.loc(:,1),d.loc(:,2),1000*ey(ind,:)',1000*ex(ind,:)','color','r','reference',refScale,'arrowstyle','tail','linewidth',2);

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled');
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(char(b(1).times(ind),'yyyy-MMM-dd HH:mm:ss'));
   
[~,indsort] = sort(abs(GIC_Subs(:,ind)),'ascend');
sortgic = GIC_Subs(indsort,ind);
sortgic(sortgic==0) = NaN;

subplot(1,2,2)
scatterm(subLoc(indsort,1),subLoc(indsort,2),3*abs(sortgic)*30,3*sortgic','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
%textm(subLoc(indsort2(1:10),1),subLoc(indsort2(1:10),2),num2cell(1:10))

bmk = [30 15 5];
scatterm([59 58 57],[-119 -119 -119],bmk*30,bmk','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'Total Peak GIC (A)';
textm([59 58 57],[-118.5 -118.5 -118.5],num2cell(bmk))
%h = scatterm(subLoc(:,1),subLoc(:,2),3*abs(GIC_Subs(:,indmaxModel)*10)+10,3*(GIC_Subs(:,indmaxModel)'),'filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%scatterm(subLoc(indnotzero,1),subLoc(indnotzero,2),abs(TneutGIC(indnotzero,indpeak))*4,TneutGIC(indnotzero,indpeak)','filled','MarkerEdgeColor','k'); hcb =colorbar; hcb.Label.String = 'GIC (A)';
%title(['Peak Neutral GIC'])
caxis([-30 50])
colormap(parula(8))
%caxis([-0.1 0.1])
set(gca,'FontSize',16)


%%
ebearing = atan2(ey,ex)*180/pi;

%plot(b(1).times,1000*max(mage(:,1:end-20),[],2),'-','Color',[0.5 0.5 0.5],'LineWidth',0.5); hold on
%plot(b(1).times,1000*min(mage(:,1:end-20),[],2),'-','Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot(b(1).times,median(1000*mage(:,1:end-20),2),'-k','LineWidth',1); hold on;
%plot(b(1).times,avg_strike,'-k','LineWidth',1); hold on;

%%
%set(gca,'YScale','log')

edges = 10.^[-7:0.1:3];
edges = [-180:5:180];

N = nan(b(1).nt,length(edges)-1);
for i = 1:b(1).nt
    %N(i,:) = histcounts(1000*mage(i,:),edges);
    N(i,:) = histcounts(ebearing(i,:),edges);
    i
end

avg_strike = nan(b(1).nt,1); std = avg_strike;
for i = 1:b(1).nt
    [avg_strike(i),std(i)] = strike_stats(ebearing(i,:),2);
    i
end

%%
pcolor(1:b(1).nt,edges(1:end-1),log10(N)'); shading flat
%%

aeso = shaperead('AESO_Planning_Areas.shp');

region = false(d.ns,4);

for i = 1:length(aeso)
    [in, on] = inpolygon(d.loc(:,2),d.loc(:,1),aeso(i).X,aeso(i).Y);
    if strcmp(aeso(i).REGION,'South') || strcmp(aeso(i).REGION,'Calgary')
        region(in==1,1) = true;
    elseif strcmp(aeso(i).REGION,'Central') || strcmp(aeso(i).REGION,'Edmonton')
        region(in==1,2) = true;
    elseif strcmp(aeso(i).REGION,'Northwest')
        region(in==1,3) = true;
    elseif strcmp(aeso(i).REGION,'Northeast')
        region(in==1,4) = true;
    end

    if strcmp(aeso(i).NAME,'Fort Sask.')
        region(in==1,4) = false;
        region(in==1,2) = true;
    end
        
end

for i = 1:4
    subplot(5,1,i+1);
    plot(b(1).times,median(1000*mage(:,region(:,i)),2),'-','LineWidth',1);
end
%%
maxe930 = max(mage(120000:130000,:),[],1);
plot(d.loc(:,2),d.loc(:,1),'.k'); hold on
indmax = find(maxe930*1000>3);
plot(d.loc(indmax,2),d.loc(indmax,1),'.r');





























%% Plots of time and magnitude of peak geoelectric field at each site
[maxe,indmaxe] = max(mage,[],1);
[sortmax,indsort2] = sort(maxe,'descend');

groups = dbscan([indmaxe; zeros(1,length(indmaxe))]',3600/2,5);

maxe = 1000*maxe;
tz = 0;
fig = figure(17);
screensize=get(groot,'Screensize');
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) 0.65*screensize(4)])

%scatter(b(1).times(indmaxe)-tz,maxe,150,maxe,'filled','MarkerEdgeColor','k'); hold on
plot(b(1).times(indmaxe)-tz,maxe,'.k'); hold on
plot(b(1).times(indmaxe(groups==1))-tz,maxe(groups==1),'ob')
plot(b(1).times(indmaxe(groups==2))-tz,maxe(groups==2),'or')
plot(b(1).times(indmaxe(groups==3))-tz,maxe(groups==3),'og')
plot(b(1).times(indmaxe(groups==4))-tz,maxe(groups==4),'oy')
%text(GIC_Times(indmaxSub(indsort2([1:10])))-tz,othermax(indsort2([1:10])),num2cell([1:10]));
text(b(1).times(indmaxe(indsort2([1:10])))-tz+hours(0.5),maxe(indsort2([1:10])),{d.site{indsort2(1:10)}});
caxis([50 60]); 
colormap(gray(5))
hcb =colorbar('Location','southoutside'); hcb.Label.String = 'Substation Latitude';
%xlim(datelimits-tz); 
grid on
set(gca,'FontSize',18)
set(gca,'layer','bottom')
ylabel('Total Peak GIC (A)')
box on


plot_lim = [47.5 61 -120.5 -109];
provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);
plot_zone = true;


%Figure 1a: 1-D Impedance
%%
fig=figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
setm(gca,'FontSize',18)

geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)
geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0,'LineWidth',2)

plotm(d.loc(groups==1,1),d.loc(groups==1,2),'.b','MarkerSize',14)
plotm(d.loc(groups==2,1),d.loc(groups==2,2),'.r','MarkerSize',14)
plotm(d.loc(groups==3,1),d.loc(groups==3,2),'.g','MarkerSize',14)
plotm(d.loc(groups==4,1),d.loc(groups==4,2),'.m','MarkerSize',14)
plotm(d.loc(groups>4,1),d.loc(groups>4,2),'.w','MarkerSize',14)

plotm([b.lat],[b.lon],'ok','MarkerSize',12,'MarkerFaceColor','w')
textm([b.lat],[b.lon],{b.site})

siteind = 1:d.ns;
indrge = min(indmaxe(groups==1)):max(indmaxe(groups==1));

[~,ind] = max(mage(indrge,:),[],1);
ind = indrge(ind);
ind = median(indrge);
ind = find(b(1).times==datetime(2024,10,10,15,15,45));
refScale = 1/10;

setm(gca,'FontSize',16)
%quivermc(d.loc(siteind,1),d.loc(siteind,2),diag(1000*ey(ind,siteind)),diag(1000*ex(ind,siteind)),'color','r','reference',refScale,'arrowstyle','tail','linewidth',2);
quivermc(d.loc(siteind,1),d.loc(siteind,2),(1000*ey(ind(1),siteind))',(1000*ex(ind(1),siteind))','color','k','reference',refScale,'arrowstyle','tail','linewidth',2);

quivermc(60.5,-118.5,5*refScale,0,'color','k','reference',refScale,'arrowstyle','tail','linewidth',3,'linestyle','filled','maxheadsize',0.01);
textm(60.8,-117.5,[num2str(5*refScale),' V/km'],'HorizontalAlignment','center')
title(char(b(1).times(ind(1))-tz,'yyyy-MMM-dd HH:mm:ss'));

gs = scaleruler;
setm(gs,'FontSize',16,'MajorTick',[0 100])





