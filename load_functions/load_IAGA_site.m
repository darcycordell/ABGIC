function b = load_IAGA_site(sites,magfile)
% Input is name of IAGA format geomagnetic field data
% Output:
%   B fields are in nT (Bx,By,Bz)
%   times: output times are a datetime object
%   site: Site name
%   lat, lon: Latitude and longitude
%
%
% NRCAN Data:

% USGS Data (API): 
%       Sitka: https://geomag.usgs.gov/ws/data/?id=SIT&starttime=2017-09-08T00:00:00Z&endtime=2017-09-08T23:59:59Z&sampling_period=1&type=variation
%       Newport: https://geomag.usgs.gov/ws/data/?id=NEW&starttime=2017-09-08T00:00:00Z&endtime=2017-09-08T23:59:59Z&sampling_period=1&type=variation

%%

Bx_all = []; By_all = []; Bz_all = []; times_all = [];
for i = 1:length(magfile)

    if strcmp(magfile{i}(1:3),sites)
        fid = fopen(magfile{i},'r');

        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        t = strsplit(fgetl(fid));
        site = t{4}; %Site name
        t = strsplit(fgetl(fid));
        lat = str2num(t{4}); %Site latitude
        t = strsplit(fgetl(fid));
        lon = str2num(t{4}); %Site longitude
        if lon>180
            lon = lon-360; %Note: IAGA longitude is 0>x>360 rather than -180>x>180
        end

        %Load Data
        while 1   
            t = fgetl(fid);
            if strcmp(t(1:4),'DATE')
                break
            end
        end
        C = textscan(fid, '%s %s %d %f %f %f %f');
        fclose(fid);

        %Load B fields
        Bx = cell2mat(C(4)); % Ignore conversion for now
        By = cell2mat(C(5));
        Bz = cell2mat(C(6));
        
        %Very slow
        %timetable = cell2table(C{1});
        %times_all = datetime(timetable.Var1,'InputFormat','yyyyMMddHHmmss');
        
        %Faster but assumes no gaps...
        t1 = datetime([C{1}{1},'-',C{2}{1}],'InputFormat','yyyy-MM-dd-HH:mm:ss.sss');
        t2 = datetime([C{1}{end},'-',C{2}{end}],'InputFormat','yyyy-MM-dd-HH:mm:ss.SSS');
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


%Remove "bad" spikes where the derivative is >0.05
indx = find((Bx_all(1:end-1)-Bx_all(2:end))./Bx_all(1:end-1)>0.05);
indy = find((By_all(1:end-1)-By_all(2:end))./By_all(1:end-1)>0.05);
ind = unique([indx; indy]);


%Remove bad points and interpolate
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