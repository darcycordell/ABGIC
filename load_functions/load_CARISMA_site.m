function [b] = load_CARISMA_site(sites,magfile)
% Input is name of CARISMA format geomagnetic field data
% Output:
%   B fields are in nT (Bx,By,Bz)
%   times: output times are a datetime object
%   site: Site name
%   lat, lon: Latitude and longitude
%
%
% CARISMA Data:
%       http://data.carisma.ca/FGM/1Hz/
%%

Bx_all = []; By_all = []; Bz_all = []; times_all = [];
for i = 1:length(magfile)

    if strcmp(magfile{i}(9:12),sites)
        fid = fopen(magfile{i},'r');

        t = strsplit(fgetl(fid));
        
        lat = str2double(t{2}); %Site latitude
        lon = str2double(t{3}); %Site longitude
        if lon>180
            lon = lon-360; %Note: IAGA longitude is 0>x>360 rather than -180>x>180
        end

        %Load Data
        C = textscan(fid, '%s %f %f %f %s');
        fclose(fid);

        %Load B fields
        Bx = cell2mat(C(2)); % Ignore conversion for now
        By = cell2mat(C(3));
        Bz = cell2mat(C(4));

        %Very slow
        %timetable = cell2table(C{1});
        %times_all = datetime(timetable.Var1,'InputFormat','yyyyMMddHHmmss');
        
        %Faster but assumes no gaps...
        t1 = datetime(C{1}(1),'InputFormat','yyyyMMddHHmmss');
        t2 = datetime(C{1}(end),'InputFormat','yyyyMMddHHmmss');
        times = t1:seconds:t2;

        N = 1:length(Bx);
        
        %Remove missing data. IAGA data seem to have missing data set with
        %99999. 
        Bx(Bx==99999)=NaN; %Set missing data to NaN
        Bx = resample(Bx,N); %Remove missing data using linear interpolation
        By(By==99999)=NaN;
        By = resample(By,N);
        Bz(Bz==99999)=NaN;
        Bz = resample(Bz,N);
        
        Bx_all = [Bx_all; Bx];
        By_all = [By_all; By];
        Bz_all = [Bz_all; Bz];
        times_all = [times_all times];
        
    end

end

[times_all, indsort] = sort(times_all);

%Set b structure
b.x = Bx_all(indsort);
b.y = By_all(indsort);
b.z = Bz_all(indsort);
b.times = times_all;
b.lat = lat;
b.lon = lon;
b.site = upper(sites);
b.nt = length(b.x);