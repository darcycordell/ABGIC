function [Zconj, Zq] = interpolate_Z(Z,freq,f,fAxis)
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

if mod(length(fAxis),2)==0
    Zconj = zeros(length(f)*2,4)+1i*zeros(length(f)*2,4);
else
    Zconj = zeros(length(f)*2-1,4)+1i*zeros(length(f)*2-1,4);
end

for ir = 1:4 %Loop through all 4 components
    % need to interpolate Z(f) onto the same set of frequencies as the B-field
    % need to do real and imaginary parts separately

    %Option #1: Interpolate onto real and imaginary (original way of doing
    %it)
%     Zqr = interp1(freq,real(Z(:,ir)),abs(fAxis(fAxis<0)),'linear')';
%     Zqi = interp1(freq,imag(Z(:,ir)),abs(fAxis(fAxis<0)),'linear')';
%     Zq(:,ir) = Zqr + 1i.*Zqi;
    
    % make values at negative frequencies the complex conjugate
    % identical results even if Z(Nyq) = 0 
    %Zconj(:,ir) = cat(1,Zq(1:end,ir),flipud(conj(Zq(1:end,ir))));
    %Zconj(:,ir) = cat(1,0,Zq(1:end,ir),flipud(conj(Zq(2:end,ir))));

    %Option #2: Better to interpolate the magnitudes and phases so
    %that you can interpolate logarithmically on mags and linearly on
    %phases. Much less error in this case.
    % This also has the benefit of resulting in a perfect extrapolation for a halfspace.
    % Option to extrapolate as: ... 'linear','extrap');
    %       However, serious danger when extrapolating to long periods
    %       because you can get some really wild extrapolations if your
    %       longer periods are a bit noisy.
    Zqmag = interp1(log10(freq),log10(abs(Z(:,ir))),log10(abs(fAxis(fAxis<0))),'linear')';
    
    logrho = interp1(log10(freq),log10((1./(2*pi*freq*4*pi*10^-7)).*abs(Z(:,ir)).^2),log10(abs(fAxis(fAxis<0))),'linear')';
    rho = 10.^logrho;
    Zqmag = log10(sqrt(rho.*2*pi.*abs(fAxis(fAxis<0))'*4*pi*10^-7));

    Zqpha = interp1(freq,angle(Z(:,ir)),abs(fAxis(fAxis<0)),'linear')';    
    
    Zq(:,ir) = 10.^(Zqmag).*exp(1i*Zqpha);

    if mod(length(fAxis),2)==0
        Zconj(:,ir) = fftshift(cat(1,conj(Zq(:,ir)),0,flip(Zq(2:end,ir))));   
    else
        Zconj(:,ir) = fftshift(cat(1,conj(Zq(:,ir)),0,flip(Zq(:,ir))));  
    end
end





