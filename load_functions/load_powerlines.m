function [lines,line_name] = load_powerlines(kmlfile)
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
count = 1; counter = 1;
while 1
    
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    tline = strtrim(tline);
    
    if strcmp(tline(1:6),'<name>')
        name = strsplit(tline,{'>','<'});
        line_name{counter} = name{3};
        counter = counter+1;
    end
    
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

segLength = 1000; %For 2021 GIC paper, I use 5000 m.

%Resample the transmission line segments into 5 km segments.
for i = 1:length(line)
    
    %Get individual transmission line path
    path = line{i}';

    %Get line segments of original transmission line
    step = [0; distance(path(1:end-1,1),path(1:end-1,2),path(2:end,1),path(2:end,2),referenceEllipsoid('earth'))];
    
    %Re-sample each segment into multiple segLength segments.
    lineseglength = cumsum(step);
    
    %If you only do linspace with the start and end points, then you end
    %up with an empty vector if lineseglength < segLength
    Vq = unique([0 linspace(0,lineseglength(end), floor(lineseglength(end)/segLength)) lineseglength(end)]);
    
    %Re-sample each segment into 100 equally spaced segments.
    %Vq = linspace(0,lineseglength(end),100);
    
    
    interpline = interp1(lineseglength, path, Vq);
    
    %New line with vertices spaced every segLength meters.
    lines{i} = interpline';


end

line_name(1:2) = [];
