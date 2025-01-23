function [DST,DSTdate] = load_dst

%https://wdc.kugi.kyoto-u.ac.jp/dstae/index.html
curdir = pwd;
cd('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/000_DATA/MAG_DATA');

fid = fopen('Dst_2010_2023.txt');
count = 1; DST = []; DSTdate = NaT;

tline = strsplit(fgetl(fid),' ');

while 1

    eof = fgetl(fid);
    if eof == -1
        break
    end

    tline = strsplit(eof,' ');

    while ~strcmp(tline{1},'DATE')
    
        tline = strsplit(fgetl(fid),' ');
    
    end
    
    while 1
    
        eof = fgetl(fid);

        if eof == -1
            break
        end

        tline = strsplit(eof,' ');
    
        if strcmp(tline{2},'Format')
            break
        end

        DSTdate(count) = datetime(string([tline{1} tline{2}]),'InputFormat','yyyy-MM-ddHH:mm:ss.SSS');
        DST(count) = str2double(tline{4});
    
        count = count+1;
    
    end
    count



end

cd(curdir);
