function b = load_raw_IAGA_site(magfile)
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


        fid = fopen(magfile,'r');

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

%%
%         %Get datetime cells
%         p=C{1}';
%         q=C{2}';
%         for i=1:length(q)
%             bb=cell2mat(q(1,i));
%             str=[' ',sprintf('%s',bb)];
%             p(1,i)=strcat(p(1,i),str);
%         end
        
        
        t1 = datetime([C{1}{1},' ',C{2}{1}],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        t2 = datetime([C{1}{end},' ',C{2}{end}],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        
        
        time = t1:seconds(1):t2;
        
        if length(time) ~= length(C{1})
            error('Something missing in your time vector!')
        end
        
        
        
%%
%        time = datetime(p, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

%%
        %Load B fields
        Bx = cell2mat(C(4)); % Ignore conversion for now
        By = cell2mat(C(5));
        Bz = cell2mat(C(6));

        N = 1:length(Bx);
        
        %Remove missing data. IAGA data seem to have missing data set with
        %99999. 
        Bx(Bx==99999)=NaN; %Set missing data to NaN
        Bx = resample(Bx,N); %Remove missing data using linear interpolation
        By(By==99999)=NaN;
        By = resample(By,N);
        Bz(Bz==99999)=NaN;
        Bz = resample(Bz,N);
        




% %Remove "bad" spikes where the derivative is >0.05
% indx = find((Bx_all(1:end-1)-Bx_all(2:end))./Bx_all(1:end-1)>0.05);
% indy = find((By_all(1:end-1)-By_all(2:end))./By_all(1:end-1)>0.05);
% ind = unique([indx; indy]);
% 
% %Meanook site has a large chunk of "bad" data from 15:09 to 15:14 UT which
% %is not a "spike" so is missed by the above derivative method.
% if strcmp(site,'MEA') 
%     ind = [ind; [54580:54860]'-tstart];
% end
% 
% %Remove bad points and interpolate
% if ~isempty(ind)
%     disp([upper(sites),': despiked ',num2str(length(ind)),' bad points'])
%     
%     Bx_all(ind) = NaN; 
%     By_all(ind) = NaN; 
%     Bz_all(ind) = NaN;
%     
%     %Interpolate NaNs
%     Bx_all = inpaint_nans(Bx_all);
%     By_all = inpaint_nans(By_all);
%     Bz_all = inpaint_nans(Bz_all);
%     
% end

%Set b structure
b.x = Bx;
b.y = By;
b.z = Bz;
b.times = time;
b.lat = lat;
b.lon = lon;
b.site = upper(site);
b.nt = length(b.x);