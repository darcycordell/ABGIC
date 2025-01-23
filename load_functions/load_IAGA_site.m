function b = load_IAGA_site(sites,magfile,sample_rate)
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
        if isempty(lat)
            lat = str2num(t{3}); 
        end
        t = strsplit(fgetl(fid));
        lon = str2num(t{4}); %Site longitude
        if isempty(lon)
            lon = str2num(t{3});
        end

        if lon>180
            lon = lon-360; %Note: IAGA longitude is 0>x>360 rather than -180>x>180
        end

        %Load Data
        formatFlag = 0;
        while 1   
            t = fgetl(fid);
            tstrspl = strsplit(t);
            if strcmp(tstrspl{3},'HDZF')
                formatFlag = 1;
            end

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
        
        if sample_rate == 1
            times = t1:seconds:t2;
        elseif sample_rate == 1/60
            times = t1:minutes:t2;
        elseif sample_rate ==1/3600
            times = t1:hours:t2;
        else
            error('Your magnetic data sample rate needs to be either second (1 Hz), minute (1/60 Hz), or hour (1/3600 Hz)')
        end

        N = 1:length(Bx);


        
        %Remove missing data. IAGA data seem to have missing data set with
        %99999. 

        Bx(Bx==99999)=NaN; %Set missing data to NaN
        By(By==99999)=NaN;
        Bz(Bz==99999)=NaN;


        if ~all(isnan(Bx))
            Bx = resample(Bx,N); %Remove missing data using linear interpolation
            By = resample(By,N);         
            Bz = resample(Bz,N);
        end
        
        Bx_all = [Bx_all; Bx];
        By_all = [By_all; By];
        Bz_all = [Bz_all; Bz];
        times_all = [times_all times];
        
    end

end

%Set b structure

if formatFlag
    %Convert from HDZF to XYZ if necessary
    b.x = Bx_all.*cos(By_all*pi/(60*180));
    b.y = Bx_all.*sin(By_all*pi/(60*180));
else
    b.x = Bx_all;
    b.y = By_all;
end

b.z = Bz_all;
b.times = times_all;
b.lat = lat;
b.lon = lon;
b.site = upper(sites);
b.nt = length(b.x);


