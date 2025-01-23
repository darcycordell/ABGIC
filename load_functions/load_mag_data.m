function [b] = load_mag_data
%%
curdir = pwd;

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
    magfile{i}
      
        if strcmp(magfile{i}(end-2:end),'sec') %NRCan

            if ~any(strcmpi(siteNames,magfile{i}(1:3))) %check if the site has already been loaded
    
                sample_rate = 1; %sample rate in Hz
                b(i) = load_IAGA_site(magfile{i}(1:3),magfile,sample_rate);
                
                siteNames = {b(:).site};
                
            end


        elseif strcmp(magfile{i}(end-2:end),'F01') %CARISMA
            
            if ~any(strcmpi(siteNames,magfile{i}(9:12))) %check if the site has already been loaded

                b(i) = load_CARISMA_site(magfile{i}(9:12),magfile);
                siteNames = {b(:).site};
            
            end
            
        elseif strcmp(magfile{i}(end-2:end),'min') %IAGA format per minute
            
            if ~any(strcmpi(siteNames,magfile{i}(1:3))) %check if the site has already been loaded
    
                sample_rate = 1/60; %sample rate in Hz
                b(i) = load_IAGA_site(magfile{i}(1:3),magfile,sample_rate);
                
                siteNames = {b(:).site};
                
            end

        elseif strcmp(magfile{i}(end-2:end),'MAG') %CANOPUS
            
            if ~any(strcmpi(siteNames,magfile{i}(9:12))) %check if the site has already been loaded

                b(i) = load_CANOPUS_site_MAG(magfile{i}(9:12),magfile);
                siteNames = {b(:).site};
            
            end


        else
            disp([magfile{i},' is an unknown file format'])  
        end
    
    
end
 

%ns = length(b); %Number of Mag sites loaded
toc

%Delete duplicates
inddelb = find(cellfun(@isempty,({b(:).x})));
b(inddelb) = [];


cd(curdir)