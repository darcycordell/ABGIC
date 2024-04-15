function [lines, line_name, line_lengths, voltage, resistance] = post_process_lines
%This function loads transmission lines from the Alberta Interconnected
%Electric System (AIES) transmission network based on publicly-available data from the Alberta
%Electric System Operator (AESO)
%
% See here: https://aeso.maps.arcgis.com/apps/webappviewer/index.html?id=959d842b42544cac9035839380b68413
%
% Additional data taken from Open Street Maps to produce the KML file in
% Google Earth. Data were manually processed and cleaned up in Google
% Earth, then some additional post-processing is required here as well.


[lines,line_name] = load_powerlines('AESO_240kV.kml');

%Add line 1114L to connect Bennett to Langdon
lines{end+1} = [50.957935 50.957133 50.956425; -113.716334 -113.719203 -113.719192];
line_name{end+1} = '1114L';

voltage = 240*ones(length(lines),1);

[lines500,line_name500] = load_powerlines('AESO_500kV.kml');
line_name500([4:8,17:18]) = [];
%%
%post-processing on 500 kV lines need to stitch together Line 13L50 and
%1325L which were in pieces in the original KML file

%Line 13L50 is in 5 pieces. In lines500 variable, the indices of the five
%pieces are 7, 3, 6, 5, 4 (in order from south to north). 

Line13L50 = [lines500{7} lines500{3} lines500{6} lines500{5} lines500{4}];

%Line 1325 is in 2 pieces. In lines500 variable, the indices of the 2
%pieces are 16 and 15 (from south to north)

Line1325 = [lines500{16} lines500{15}];

lines500{3} = Line13L50;
lines500{15} = Line1325;

lines500([4 5 6 7 16]) = [];

%% After post-processing, we can take 500 kV and 240 kV lines and combine them into the same cell structure

lines = [lines lines500];
line_name = [line_name line_name500];

voltage = [voltage; 500*ones(length(lines500),1)];

%% Need to sort lines so that each line starts from the north and goes to the south

for i = 1:length(lines)
    firstlat = lines{i}(1,1);
    lastlat = lines{i}(1,end);

    if firstlat < lastlat
        lines{i}(1,:) = lines{i}(1,end:-1:1);
        lines{i}(2,:) = lines{i}(2,end:-1:1);
    end


end

%% Get Line lengths
line_lengths = nan(length(lines),1); resistance = nan(length(lines),1);
for i = 1:length(lines)
    line_lengths(i) = sum(distance(lines{i}(1,1:end-1),lines{i}(2,1:end-1),lines{i}(1,2:end),lines{i}(2,2:end),referenceEllipsoid('wgs84')));
    if voltage(i) == 240
        %I am using 0.0455 Ω/mile for 240 kV lines (same as Horton et al.
        % uses for 345 kV lines).
        resistance(i) = 0.0455*line_lengths(i)/1.60934/1000; %convert from meters to miles
    elseif voltage(i) == 500
        %I am using 0.0227 Ω/mile for 500 kV lines (same as Horton et al.
        % uses for 500 kV lines).
        resistance(i) = 0.0227*line_lengths(i)/1.60934/1000; %convert from meters to miles
    else
        resistance(i) = NaN;
    end
end

