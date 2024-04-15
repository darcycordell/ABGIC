function [xt,yt] = calc_ifft(X,Y,pad)
%Calculate ifft for X and Y frequency-domain vector data
%
% Inputs:
%       X = frequency-domain x component
%       Y = frequency-domain y component
%       pad = zero padding
%
% Outputs:
%       xt = time-domain x component
%       yt = time-domain y component
%       f = vector of frequencies

if all(isnan(X(2:end))) || all(isnan(Y(2:end)))
    xt = zeros(length(X),1);
    yt = zeros(size(xt));
else

    N = 1:length(X);

    X = resample(X,N);
    Y = resample(Y,N);

    xt = ifft(X,'symmetric');
    yt = ifft(Y,'symmetric');

end

%Remove padding
xt(1:pad) = [];
xt(end-pad+1:end) = [];
yt(1:pad) = [];
yt(end-pad+1:end) = [];
