function [predX,predY,predZ,LON,LAT,JX,JY]=interpolate_secs(b,lim,dlon,dlat)
% Option to set lim = NaN and interpolate on the (dlat,dlon) points instead
tic
obsX = [b(:).x];
obsY = [b(:).y];
obsZ = [b(:).z];

obsLat = [b(:).lat];
obsLon = [b(:).lon];

baselineX = 0;%nanmedian(obsX);
baselineY = 0;%nanmedian(obsY);
baselineZ = 0;%nanmedian(obsZ);

if ~isnan(lim)
    %Query points for B-field interpolated
    [qLAT,qLON] = meshgrid(lim(1):dlat:lim(2),lim(3):dlon:lim(4));
    
    %Points where SECS are calculated. Staggered from B-field query points, and
    %with grid expanded by 2*dlat and 2*dlon
    [SECSlat,SECSlon] = meshgrid(lim(1)-dlat*2-dlat/2:dlat:lim(2)+dlat*2+dlat/2,lim(3)-dlon*2-dlon/2:dlon:lim(4)+dlon*2+dlon/2);
else
    qLAT = dlat;
    qLON = dlon;
    [SECSlat,SECSlon] = meshgrid(linspace(min(dlat)-10,max(dlat)+10,50),linspace(min(dlon)-10,max(dlon)+10,50));

end

LAT = qLAT';
LON = qLON';

[predX,predY,predZ,JX,JY] = secs(obsLat,obsLon,obsX-baselineX,obsY-baselineY,obsZ-baselineZ,SECSlat(:),SECSlon(:),LAT(:),LON(:));



toc
