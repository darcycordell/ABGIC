function [b] = load_mag_data(s)
curdir = pwd;

[magfile, magpath]=uigetfile({'*'},'Choose Mag Files','MultiSelect','on');
if ~iscell(magfile)
    if magfile == 0
        return
    end
end
cd(magpath);

t1 = datetime(2017,9,8,0,0,0);
t2 = datetime(2017,9,8,23,59,59);
times = t1:seconds(1):t2;

tic
for i = length(magfile):-1:1
    if strcmp(magfile{i}(end-2:end),'sec')
        b(i) = load_IAGA_site(magfile{i}(1:3),magfile,times);
        
    elseif strcmp(magfile{i}(end-2:end),'F01')
        b(i) = load_CARISMA_site(magfile{i}(9:12),magfile,times);
        
    else
        disp([magfile{i},' is an unknown file format'])  
    end
    
end
 

ns = length(b); %Number of Mag sites loaded
toc
%Option to include synthetic signal (sine wave with 1000 nT peak at variable frequency)
if s.flag
    %If you are doing synthetic data then you set a synthetic frequency
    for i = 1:ns
        b(i).times = b(i).times(1:s.len);
        
        %Polarized in North (x) direction
        b(i).x = s.mag*cos(2*pi*s.freq*[0:length(b(i).times)-1])';
        b(i).y = zeros(size(b(i).x));
        b(i).z = zeros(size(b(i).x));
        b(i).nt = length(b(i).x);
    end
end

cd(curdir)