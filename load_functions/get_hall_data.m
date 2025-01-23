function [GIC_Hall,t] = get_hall_data


get_files;

if strcmp(hallFile,'none')
    GIC_Hall = 0;
    t = 0;
    return
end

curdir = pwd;
cd(hallFilePath)
%GIC_Hall = readtable('AltaLink_GIC_230322_25_merge.xlsx');
GIC_Hall = readtable(hallFile);
%GIC_Hall = readtable('AltaLink_GIC_211012.csv');
cd(curdir)

%The raw Excel data does not have seconds
tgic = datetime(GIC_Hall{:,1}+hours(6)-seconds(0), 'Format','yyyy-MM-dd HH:mm:ss');
%%

if ~any(second(tgic)~=0)
%Find first complete minute
t1 = find(tgic==tgic(1)+minutes(1));
if length(t1)>30
    t1 = [];
end
%Find last complete minute
t2 = find(tgic==tgic(end)-minutes(1));

%First indices of leading seconds
beginning = find(tgic==tgic(1));

%Last indices of trailing seconds
ending = find(tgic==tgic(end));

firstTime = tgic(t1(1))-seconds(2*length(beginning));
lastTime = tgic(t2(end))+minutes(1)+seconds(2*(length(ending)-1));

%Hall probe is sampled every 2 seconds from the first time to the last time
t = (firstTime:seconds(2):lastTime)';


if length(t) ~= length(tgic)
    t = (firstTime:seconds(2):lastTime-1)';
    warning('Your Hall probe data does not have equally spaced data!!')
end

else
    t = tgic; %data includes seconds already
end
    



