function [fwd]=calc_fwd_1d(model_depth,model_res,freq_array,err)
%
% Function which calculates the 1D MT data from a given 1D model using
% Waits recursion. The model should be in the following format:
%
% Depth (m)         Resistivity (Ohm m)
%  depth 1              resistivity 1
%  depth 2              resistivity 2
%   ...                     ...
%
% Depth is defined as the depth to the top of the layer!
%
% Usage: [fwd] = calc_fwd_1d(model_depth,model_res,freq_array,err)
%
% Inputs:
%
%   model_depth: vector of depths (to top of layer)
%   model_res: vector of resistivities of each layer
%   freq_array: vector of frequencies
%   err: Error level (adds % Gaussian noise)
%
% Outputs:
%
%   fwd: A structure containing the apparent resistivity, phase, impedance
%   and admittance for the 1D model.

%Set up constants given a frequency array
mu=4*pi*10^-7;
w=(2*pi).*freq_array; %angular frequency
sigma=1./model_res; %conductivity (S/m)
num=length(model_res);
nf = length(freq_array);

%Calculate the wave number and admittance values for the halfspace
k(num,:)=sqrt(-1i.*mu.*w.*sigma(num));
C(num,:)=1./k(num,:);

%Wait's recursion relation to calculate wave number and admittance of each layer
for n=(num-1):-1:1
    k(n,:)=sqrt(-1i.*w.*mu.*sigma(n));
    C(n,:)=(1./k(n,:)).*(C(n+1,:).*k(n,:)+tanh(k(n,:).*(model_depth(n+1)-model_depth(n))))./ ... 
        (C(n+1,:).*k(n,:).*tanh(k(n,:).*(model_depth(n+1)-model_depth(n)))+1);
end

%Impedance, apparent resistivity and phase
Zxy=w.*mu.*-conj(C(1,:))*exp(1i*-pi/2);

C = conj(C(1,:));

if isnan(err)
    %if error isnan then use "variable error" as per Equation 4 of the
    %Supplementary Material of GRL Manuscript
    a = 4-2*log(0.19);
    err = exp(0.5*(log10(1./freq_array)-a))+0.01;
end

Zer = abs(Zxy).*randn(1,nf).*err;
Zei = abs(Zxy).*randn(1,nf).*err;

rZ = real(Zxy)+Zer;
iZ = imag(Zxy)+Zei;

Zxy = rZ + 1i.*iZ; %Impedance


fwd.rho=(1./(w.*mu)).*(abs(Zxy)).^2;
fwd.phi=(180./pi).*atan2(imag(Zxy),real(Zxy));
fwd.Z = Zxy;
fwd.Zerr = abs(Zxy).*err;
fwd.C = C;
fwd.freq_array = freq_array;

end


