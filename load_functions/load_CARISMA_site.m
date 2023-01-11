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
        
        Bx_all = [Bx_all; Bx];
        By_all = [By_all; By];
        Bz_all = [Bz_all; Bz];
        times_all = [times_all times];
        
    end

end


%Remove bad points and interpolate
ind = unique([find(Bx_all>=99999.9); find(By_all>=99999.9); find(Bz_all>=99999.9)]);
if ~isempty(ind)
    disp([upper(sites),': despiked ',num2str(length(ind)),' bad points'])

    Bx_all(ind) = NaN;
    By_all(ind) = NaN;
    Bz_all(ind) = NaN;
    
    %Interpolate NaNs
    Bx_all = inpaint_nans(Bx_all);
    By_all = inpaint_nans(By_all);
    Bz_all = inpaint_nans(Bz_all);
    
end

%Set b structure
b.x = Bx_all;
b.y = By_all;
b.z = Bz_all;
b.times = times_all;
b.lat = lat;
b.lon = lon;
b.site = upper(sites);
b.nt = length(b.x);