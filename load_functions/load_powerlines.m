function lines = load_powerlines(kmlfile)
%
% Function which loads power transmission lines from a Google Earth KML
% file. The file (kmlfile) must be on path or in current directory.
% The raw lines are then cut into 5 km segments
%
% Outputs: lines is a cell array whose entries are transmission lines
%       Each entry is a (2 x nseg) matrix of lat-long pairs for a given
%       transmission line. The lines are resampled into consistent 5 km
%       segments.
%
%
fid = fopen(kmlfile);

nlon = 0;
nlat = 0;
count = 1;
while 1
    
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    tline = strtrim(tline);
    
    if strcmp(tline,'<coordinates>')
        tline = fgetl(fid);
        dat = strsplit(tline,{' ',','});
        dat = str2double(reshape(dat(1:end-1),3,(length(dat)-1)/3));
        
        %Raw lat-long vertices of the transmission lines
        % The raw transmission lines usually only have 2 - 10 vertices
        % which are 50 - 100 km long each. This is not finely sampled
        % enough for the geoelectric field line integrals.
        line{count} = dat([2 1],:);
        
        count = count+1;
    end
        
    
end

fclose(fid);

%Resample the transmission line segments into 5 km segments.
for i = 1:length(line)
    
    %Get individual transmission line path
    path = line{i}';

    %Get line segments of original transmission line
    step = [0; distance(path(1:end-1,1),path(1:end-1,2),path(2:end,1),path(2:end,2),referenceEllipsoid('earth'))];
    
    %Re-sample each segment into multiple 5 km segments.
    lineseglength = cumsum(step);
    Vq = linspace(0,lineseglength(end), round(lineseglength(end)/5000));
    interpline = interp1(lineseglength, path, Vq);
    
    %New line with vertices spaced very 5 km.
    lines{i} = interpline';


end
