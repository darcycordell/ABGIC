function [b] = load_CANOPUS_site_MAG(sites,magfile)
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


        tmp = fgetl(fid);
        while strcmp(tmp(1),'#')
            tmp= fgetl(fid);
        end

%         if fid==-1
%             b = [];
%             return
%         end
%         line = fgetl(fid);
% 
%         if line == -1
%             b = [];
%             return
%         end

        t = strsplit(tmp);
        
        lat = str2double(t{2}); %Site latitude
        lon = str2double(t{3}); %Site longitude
        if lon>180
            lon = lon-360; %Note: IAGA longitude is 0>x>360 rather than -180>x>180
        end

        %Load Data
        C = textscan(fid, '%s %f %f %f %s');
        fclose(fid);

        if isempty(C{1})
            b.x = nan(17280,1);
            b.y = nan(17280,1);
            b.z = nan(17280,1);
            b.times = NaT(17280,1);
            b.nt = 17280;
            disp('EMPTY FILE')
            return
        end

        %Load B fields
        Bx = cell2mat(C(2)); % Ignore conversion for now
        By = cell2mat(C(3));
        Bz = cell2mat(C(4));

        badflag = find(~strcmp(C{5},'.'));

        Bx(badflag) = NaN;
        By(badflag) = NaN;
        Bz(badflag) = NaN;

        %Remove missing data. IAGA data seem to have missing data set with
        %99999. 
        bad9 = abs(Bx-99999)<1 | abs(Bx+80692)<1;
        Bx(bad9)=NaN; %Set missing data to NaN 
        By(bad9)=NaN;
        Bz(bad9)=NaN;



        %Very slow
        %timetable = cell2table(C{1});
        %times_all = datetime(timetable.Var1,'InputFormat','yyyyMMddHHmmss');
        
        %Faster but assumes no gaps...
        t1 = datetime(C{1}{1}(1:8),'InputFormat','yyyyMMdd');
        times = t1:seconds(5):datetime(year(t1),month(t1),day(t1),23,59,59);

        if length(times) == length(Bx)+1
            if abs(times(end-1)-datetime(C{1}(end),'InputFormat','yyyyMMddHHmmss'))<5
                Bx(end+1) = NaN;
                By(end+1) = NaN;
                Bz(end+1) = NaN;
            else
                flag = 'bad';
            end
        end

        if length(Bx)~=17280
            tcheck = datetime(C{1},'InputFormat','yyyyMMddHHmmss');

            [tcheck,ind] = unique(tcheck);

            Bx = Bx(ind);
            By = By(ind);
            Bz = Bz(ind);

            tcheck = tcheck - seconds(rem(second(tcheck),5));

            [tcheck,ind] = unique(tcheck);

            Bx = Bx(ind);
            By = By(ind);
            Bz = Bz(ind);

            [~,ia,ib] = intersect(tcheck,times);
            x = nan(length(times),1);
            y = nan(length(times),1);
            z = nan(length(times),1);

            x(ib) = Bx(ia);
            y(ib) = By(ia);
            z(ib) = Bz(ia);

            Bx = x;
            By = y;
            Bz = z;

            %plot(times,Bx)
            %xlim([times(1) times(end)])

        end

        
        %Bx = resample(Bx,N); %Remove missing data using linear interpolation
        %By = resample(By,N);
        %Bz = resample(Bz,N);
        
        Bx_all = [Bx_all; Bx];
        By_all = [By_all; By];
        Bz_all = [Bz_all; Bz];
        times_all = [times_all times];


        nbad9 = sum(bad9);
        nbadflag = length(badflag);
        
    end

end

%disp(['Bad 99999 removed = ',num2str(nbad9)])
%disp(['Bad Flag removed = ',num2str(nbadflag)])

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