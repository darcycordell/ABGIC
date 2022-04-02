function Z = get_Z_at_location(lat,long,method,save_file)

% method = nearest (default), linear, or natural

if ~exist('method','var')
    method = 'nearest';
end

[zn] = load_trich_zones;
[d,in,~] = load_assign_impedance(zn,'AB_BC_MT_DATA_512_sites.mat');
d.Z(abs(real(d.Z(:)))>10^5)=NaN; %remove bad outliers

[a,~,~] = calc_Z_trich(d,zn,in);

[dist, indlatlon] = sort(distance(lat,long,d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84')));  
disp(['Interpolation Method: ',method,'. Nearest MT Site: ',d.site{indlatlon(1)},' located ', num2str(round(dist(1)/1000,1)), ' km away'])

ind = indlatlon(1:end);

Z = nan(d.nf,4);
for i = 1:d.nf
    
        Fxx =  scatteredInterpolant(d.loc(ind,2),d.loc(ind,1),squeeze(d.Z(i,1,ind)),method);
        Fxy =  scatteredInterpolant(d.loc(ind,2),d.loc(ind,1),squeeze(d.Z(i,2,ind)),method);
        Fyx =  scatteredInterpolant(d.loc(ind,2),d.loc(ind,1),squeeze(d.Z(i,3,ind)),method);
        Fyy =  scatteredInterpolant(d.loc(ind,2),d.loc(ind,1),squeeze(d.Z(i,4,ind)),method);
        
        Vxx = Fxx(long,lat);
        Vxy = Fxy(long,lat);
        Vyx = Fyx(long,lat);
        Vyy = Fyy(long,lat);
        
        Z(i,:) = [Vxx Vxy Vyx Vyy];
        
        %disp(['E-field for MT Site #: ',num2str(i),' (Name: ',d.site{i},') completed'])
end

%d = edit_data(d);
%%

[dp] = dplus(d.Z(:,2,indlatlon(1)),d.Zerr(:,2,indlatlon(1)),d.T,10);

[stats] = statistics([real(dp.Zorig); imag(dp.Zorig)],[real(dp.Z); imag(dp.Z)],[real(dp.Zerr); imag(dp.Zerr)]);

indres = find(abs(stats.r)>2.5*stats.rms);

indres(indres>length(dp.T)) = indres(indres>length(dp.T)) - length(dp.T);

subplot(1,2,1)
loglog(d.T,real(d.Z(:,2,indlatlon(1))),'xr'); hold on; grid on
title('Real XY')
%loglog(dp.T(indres),real(dp.Zorig(indres)),'xb'); hold on


[dp] = dplus(d.Z(:,3,indlatlon(1))*exp(1i*pi),d.Zerr(:,3,indlatlon(1)),d.T,10);

[stats] = statistics([real(dp.Zorig); imag(dp.Zorig)],[real(dp.Z); imag(dp.Z)],[real(dp.Zerr); imag(dp.Zerr)]);

indresyx = find(abs(stats.r)>2.5*stats.rms);

indresyx(indresyx>length(dp.T)) = indresyx(indresyx>length(dp.T)) - length(dp.T);

subplot(1,2,2)
loglog(d.T,-real(d.Z(:,3,indlatlon(1))),'xb'); hold on; grid on
title('Real YX')
%loglog(dp.T(indresyx),real(dp.Zorig(indresyx)),'xb'); hold on

%loglog(d.T,real(a.Z(:,2,indlatlon(1))),'xr');

indres = unique(indres);

%%  
% provinces = shaperead('province.shp','UseGeoCoords',true);
% states = shaperead('usastatehi','UseGeoCoords',true);
% 
% plot_lim = [min(d.loc(indlatlon(1:50),1)) max(d.loc(indlatlon(1:50),1)) min(d.loc(indlatlon(1:50),2)) max(d.loc(indlatlon(1:50),2))];
% screensize=get(groot,'Screensize');
% fig=figure(1);
% clf
% set(fig,'Position',[0.1*screensize(3) 0*screensize(4) 0.35*screensize(3) screensize(4)])
% worldmap([plot_lim(1) plot_lim(2)],[plot_lim(3) plot_lim(4)]);
% 
% geoshow(states,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0)
% geoshow(provinces,'DefaultFaceColor',[0 0 0],'DefaultEdgeColor','black','facealpha',0)
% 
% plotm(d.loc(:,1),d.loc(:,2),'.k','MarkerSize',10); hold on
% plotm(d.loc(ind,1),d.loc(ind,2),'.k','MarkerSize',18); hold on
% plotm(d.loc(indlatlon(1),1),d.loc(indlatlon(1),2),'.r','MarkerSize',15); hold on
% plotm(lat,long,'xr','MarkerSize',14,'LineWidth',2);

%%
if save_file
data_matrix = d.f;

fid = fopen([num2str(lat),'_',num2str(long),'_ImpedanceSite_',d.site{indlatlon(1)},'.txt'],'w');
type = '%-15s';
fprintf(fid,[type type type type type type type type type '\n'],'Freq','real(Zxx)','imag(Zxx)','real(Zxy)','imag(Zxy)','real(Zyx)','imag(Zyx)','real(Zyy)','imag(Zyy)');


for i = 1:4
    
    data_matrix = [data_matrix real(d.Z(:,i,indlatlon(1))) imag(d.Z(:,i,indlatlon(1)))];
    
end

data_matrix(isnan(d.Z(:,i,indlatlon(1))),:) = [];
type = '%-15.5e';
fprintf(fid,[type type type type type type type type type '\n'],data_matrix');

fclose(fid);

save([num2str(lat),'_',num2str(long),'_ImpedanceSite_',d.site{indlatlon(1)},'.mat'],'Z')
end