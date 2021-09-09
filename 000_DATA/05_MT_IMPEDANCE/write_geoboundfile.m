provinces = shaperead('province.shp','UseGeoCoords',true);
states = shaperead('usastatehi','UseGeoCoords',true);

fid = fopen('provinces.txt','w');
for i = 1:length(provinces)
    Lon = provinces(i).Lon(~isnan(provinces(i).Lon));
    Lat = provinces(i).Lat(~isnan(provinces(i).Lon));
    fprintf(fid,'\t%.5f\t%.5f\n',[Lon; Lat]);
    fprintf(fid,'%s\n','X');
end

fclose(fid);


fid = fopen('states.txt','w');
for i = 1:length(states)
    Lon = states(i).Lon(~isnan(states(i).Lon));
    Lat = states(i).Lat(~isnan(states(i).Lon));
    fprintf(fid,'\t%.5f\t%.5f\n',[Lon; Lat]);
    fprintf(fid,'%s\n','X');
end

fclose(fid);