function [bx,by,Bx,By,b] = get_b(latq,lonq)
curdir = pwd; 

get_files;


disp('*********************LOAD OBSERVARTORY DATA************************')
cd(bFilePath)
load(bFile);
cd(curdir)
disp(['Done loading magnetic observatory data from file: ',bFile])

b(1).fs = 1;
b(1).pad = 10000;


disp('...converting observatory data to frequency domain')
ns = length(b);
%Convert mag sites to the frequency domain.
b(ns).X = []; b(ns).Y = []; b(ns).f = [];
for i = 1:ns
    [b(i).X,b(i).Y,b(i).f,b(i).fAxis] = calc_fft(b(i).x,b(i).y,b(1).fs,b(1).pad);

end

nf = length(b(1).fAxis);

disp('Interpolating b(t) Spatially...')

[bx,by] = interpolate_b(b,latq,lonq);
disp('...done interpolation')
disp('...converting interpolated b to frequency domain ...')

ttic = tic;
Bx = zeros(length(latq),nf); Bx = Bx+1i*Bx;
By = Bx;
for i = 1:size(bx,1)
    [Bx(i,:),By(i,:)] = calc_fft(bx(i,:)',by(i,:)',b(1).fs,b(1).pad); 
end


disp(['...completed B(omega) in ',num2str(toc(ttic)),' seconds'])

