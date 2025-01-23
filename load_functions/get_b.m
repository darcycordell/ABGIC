function [b] = get_b(flag)
curdir = pwd; 

get_files;


disp('*********************LOAD OBSERVARTORY DATA************************')

if strcmp(flag,'mat')
    cd(bFilePath)
    load(bFile);
    cd(curdir)
elseif strcmp(flag,'raw')
    b = load_mag_data;
end

disp(['Done loading magnetic observatory data from file: ',bFile])

freq_menu = menu('Sample rate:','1 Hz','1 min')

if freq_menu == 2
    b(1).fs = 1/60;
else
    b(1).fs = 1;
end

b(1).pad = 10000;


disp('...converting observatory data to frequency domain')
ns = length(b);
%Convert mag sites to the frequency domain.
b(ns).X = []; b(ns).Y = []; b(ns).f = [];
for i = 1:ns
    [b(i).X,b(i).Y,b(i).f,b(i).fAxis] = calc_fft(b(i).x,b(i).y,b(1).fs,b(1).pad);

end




