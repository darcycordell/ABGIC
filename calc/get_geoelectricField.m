function [ex,ey] = get_geoelectricField(magfile,latq,lonq,fs,plot_flag)
% Function which computes the geoelectric field at a given magnetometer
% location. Interpolates gaps using an external function called
% inpaint_nans which solves a linear system to interpolate/extrapolate gaps
% prior to performing the fft
%
% Inputs:
%       magfile: File name for magnetic data in nT
%       latq: Latitude of the magnetometer station location
%       lonq: Longitude (+/-180°) of the magnetometer station location
%       fs = Sampling rate in Hz
%       plot_flag = true false flag to plot output
%
% Outputs:
%       ex = North-component of electric field in V/m
%       ey = East-component of electric field in V/m
%
% Examples:
%   magfile = '20141017FCHP.F01';
%   latq = 58.769; %query latitude to interpolate onto
%   lonq = -111.106; %query longitude to interpolate onto
%   fs = 1;
%
%   magfile = '19990101MCMU.MAG';
%   latq = 56.657;
%   lonq = -111.21;
%   fs = 1;

%Open and load raw mag data
fid = fopen(magfile);
C = textscan(fid, '%f %f %f %f %s');
t = datetime(datenum(num2str(C{1},'%d'), 'yyyymmddHHMMSS'),'ConvertFrom','datenum');
fclose(fid);

%Unedited bx (north) and by (east) data
bx_orig = C{2};
by_orig = C{3};
%bz = C{4};

%bad flags are x
badFlag = strcmp(C{5},'x');

%set bad points to NaN
bx_orig(badFlag) = NaN;
by_orig(badFlag) = NaN;

%interpolate over NaN gaps prior to FFT
bx = inpaint_nans(bx_orig);
by = inpaint_nans(by_orig);

nt = length(bx); %length of time series

%Load MT impedance array data (file must be on your path)
load('AB_BC_MT_DATA_526_sites_230321.mat');
d.Z(abs(real(d.Z(:)))>10^5)=NaN; %remove impedance outliers

pad = 10000; %FFT padding

bx = squeeze(bx); by = squeeze(by);

%Find the MT site index which is closest to the magnetometer site location
[~,indMT] = min(distance(latq,lonq,d.loc(:,1),d.loc(:,2),referenceEllipsoid('wgs84'))); 

nf = nt+2*pad; %Number of frequencies in FFT  including negative freqs

%Get magnetic data in frequency domain
%   f = positive frequencies only
[bX,bY,f] = calc_fft(bx,by,fs,pad);

%Frequency indices to FFT
fidx = 1:nf;
  
%Get Z impedance on the same frequency set
Zint3D =  interpolate_Z(d.Z(:,:,indMT),d.f,f);

%Compute E fields in frequency domain
[Ex3D,Ey3D] = calc_E(bX(fidx).',bY(fidx).',Zint3D(fidx,:));

%Compute E fields in time domain
[ex,ey] = calc_ifft(Ex3D,Ey3D,pad);

%Optional plot

if plot_flag
    subplot(2,1,1)
    plot(t,bx_orig-nanmean(bx_orig)); hold on;
    plot(t,by_orig-nanmean(by_orig)); grid on
    xlim([t(1) t(end)])
    legend('b_x(t)','b_y(t)')
    ylabel('Magnetic Field (nT)')
    
    subplot(2,1,2)
    plot(t,1000*ex); hold on;
    plot(t,1000*ey); grid on
    xlim([t(1) t(end)])
    legend('e_x(t)','e_y(t)')
    ylabel('Electric Field (V/km)')
end
