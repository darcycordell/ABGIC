function [b] = load_mag_data(varargin)
curdir = pwd;

if isempty(varargin)
    s.flag = false;
elseif length(varargin)==1
    s = varargin{1};
else
    error('load_mag_data: Too many input argument')
end


[magfile, magpath]=uigetfile({'*'},'Choose Mag Files','MultiSelect','on');
if ~iscell(magfile)
    if magfile == 0
        return
    end
end
cd(magpath);

if ischar(magfile)
    magfile = {magfile};
end

tic
siteNames = {''};
for i = length(magfile):-1:1
      
        if strcmp(magfile{i}(end-2:end),'sec')

            if ~any(strcmpi(siteNames,magfile{i}(1:3))) %check if the site has already been loaded
    
                sample_rate = 1; %sample rate in Hz
                b(i) = load_IAGA_site(magfile{i}(1:3),magfile,sample_rate);
                
                siteNames = [{b(:).site}];
                
            end


        elseif strcmp(magfile{i}(end-2:end),'F01')
            
            if ~any(strcmpi(siteNames,magfile{i}(9:12))) %check if the site has already been loaded

                b(i) = load_CARISMA_site(magfile{i}(9:12),magfile);
                siteNames = [{b(:).site}];
            
            end
            
        elseif strcmp(magfile{i}(end-2:end),'min') %IAGA format per minute
            
            if ~any(strcmpi(siteNames,magfile{i}(1:3))) %check if the site has already been loaded
    
                sample_rate = 1/60; %sample rate in Hz
                b(i) = load_IAGA_site(magfile{i}(1:3),magfile,sample_rate);
                
                siteNames = [{b(:).site}];
                
            end


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