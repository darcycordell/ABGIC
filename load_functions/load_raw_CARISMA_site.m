function [b] = load_raw_CARISMA_site(magfile)
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


        fid = fopen(magfile,'r');

        t = strsplit(fgetl(fid));
        
        site = t{1};
        lat = str2double(t{2}); %Site latitude
        lon = str2double(t{3}); %Site longitude
        if lon>180
            lon = lon-360; %Note: IAGA longitude is 0>x>360 rather than -180>x>180
        end

        %Load Data
        C = textscan(fid, '%s %f %f %f %s');
        fclose(fid);

        %time = datetime(C{1},'InputFormat','yyyyMMddHHmmss');
        %Get datetime cells
        
        t1 = datetime(C{1}(1),'InputFormat','yyyyMMddHHmmss');
        t2 = datetime(C{1}(end),'InputFormat','yyyyMMddHHmmss');
        
        time = t1:seconds(1):t2;
        
        if length(time) ~= length(C{1})
            error('Something missing in your time vector!')
        end
        
        %Load B fields
        Bx = cell2mat(C(2)); % Ignore conversion for now
        By = cell2mat(C(3));
        Bz = cell2mat(C(4));


        N = 1:length(Bx);
        
 
        
 



%Remove bad points and interpolate
ind = unique([find(Bx>=99999.9); find(By>=99999.9); find(Bz>=99999.9)]);


%Remove bad points and interpolate
if ~isempty(ind)
    disp([upper(site),': despiked ',num2str(length(ind)),' bad points'])
    
    Bx(ind) = NaN; 
    By(ind) = NaN; 
    Bz(ind) = NaN;
    
    %Interpolate NaNs
    Bx = inpaint_nans(Bx);
    By = inpaint_nans(By);
    Bz = inpaint_nans(Bz);
    
end

%Set b structure
b.x = Bx;
b.y = By;
b.z = Bz;
b.times = time;
b.lat = lat;
b.lon = lon;
b.site = upper(site);
b.nt = length(b.x);