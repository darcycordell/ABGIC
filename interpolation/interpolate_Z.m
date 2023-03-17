function [Z_int, Zq] = interpolate_Z(Z,freq,f)
%Interpolate Z onto the same frequency set as the magnetic B-field data
%Also do some manipulation (e.g. conjugation, units etc.) to prepare
% Z for the calculation of the frequency-domain E-field
%
% Inputs: Z (nf x 4) impedance tensor with xx, xy, yx, yy data for a single
%               site location
%         freq is the original frequency set from the MT data
%         f is the new frequency set to interpolate to from the B-field data
%
% Outputs: Z_int is (nf_new x 4) interpolated impedance
%
% Note that Z is in SI units of Ohm

indnan = isnan(Z(:,2));
Z(indnan,:) = [];
freq(indnan) = [];

Zq = zeros(length(f),4)+1i*zeros(length(f),4);
Zconj = zeros(length(f)*2,4)+1i*zeros(length(f)*2,4);
for ir = 1:4 %Loop through all 4 components
    % need to interpolate Z(f) onto the same set of frequencies as the B-field
    % need to do real and imaginary parts separately
    Zqr = interp1(freq,real(Z(:,ir)),f,'linear')';
    Zqi = interp1(freq,imag(Z(:,ir)),f,'linear')';
    Zq(:,ir) = Zqr + 1i.*Zqi;
    
    % make values at negative frequencies the complex conjugate
    % identical results even if Z(Nyq) = 0 
    %Zconj(:,ir) = cat(1,Zq(1:end,ir),flipud(conj(Zq(1:end,ir))));
    Zconj(:,ir) = cat(1,0,Zq(1:end,ir),flipud(conj(Zq(2:end,ir))));
end

Z_int = Zconj;

