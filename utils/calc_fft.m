function [X,Y,f,fAxis] = calc_fft(xt,yt,fs,pad)
%Calculate fft for x and y vector data
%
% Inputs:
%       xt = time-domain x component
%       yt = time-domain y component
%       fs = sample rate (e.g. fs = 1 is 1 Hz; fs = 1/60 is 1 minute sample
%       pad = zero padding
%
% Outputs:
%       X = frequency-domain x component
%       Y = frequency-domain y component
%       f = vector of frequencies

% add padding to beginning and end of time series
%xt = cat(1,ones(pad,1).*xt(1),xt,ones(pad,1).*xt(end));%-nanmean(xt);
%yt = cat(1,ones(pad,1).*yt(1),yt,ones(pad,1).*yt(end));%-nanmean(yt);

xt = cat(1,ones(pad,1).*xt(1),xt,ones(pad,1).*xt(end));%-nanmean(xt);
yt = cat(1,ones(pad,1).*yt(1),yt,ones(pad,1).*yt(end));%-nanmean(yt);

L = length(yt); % number of samples to fft

X = fft(xt); % fft Bx
Y = fft(yt); % fft By

% f = fs*(0:(L/2))/L;
%% 
df = fs/L;
fAxis = (0:df:(fs-df)) - (fs-mod(L,2)*df)/2;

fAxis(abs(fAxis)<df/1000000) = 0;

%Identical to fAxis for even length signal
%fAxis2 = -fs/2:df:fs/2-df;

%Different from fAxis, includes double 0 which I think is wrong
%fAxis3 = horzcat(-linspace(0,L/2,L/2)*fs/L,linspace(L/2,0,L/2)*fs/L);

%Shifted from fAxis by 1 sample
%fAxis4 = fs/L*[(0:L/2) -1*((L/2-1):-1:1)];

%Identical to fAxis for both odd and even length signals
% fAxis5 = 0 : df : df*(L-1);
% fAxis5 = fftshift(fAxis5);
% ind = fAxis5>=fs/2-eps(fs/2);
% fAxis5(ind) = fAxis5(ind)-fs;

%f = linspace(0,fs,L/2)-(fs/2);

if mod(L,2)==0
   f = fAxis((L/2)+1:end);
else
   f = fAxis((L+1)/2:end);
end

