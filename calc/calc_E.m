function [Ex,Ey] = calc_E(Bx,By,Z)
% Function to calculate the E-field vector given B-field vector and 2x2
% impedance in frequency domain
%
% Inputs: Bx and By are north (x) and east (y) magnetic fields in nT at a
%           specific location in the frequency domain
%           Size is (nf x 1)
%         Z is impedance (Ohms) at a specific location
%           Size is (nf x 4) matrix of xx, xy, yx, yy impedances at nf frequencies
%
% Outputs: Ex and Ey are north and east electric fields in V/m at the same
%           location as Bx and By

%%
mu = 4*pi*10^-7;

Z = reshape(Z(:,[1 3 2 4]).',2,2,size(Bx,2));
B = [Bx; By]*10^-9; %Conversion from nT to Tesla

B = permute(B,[1,3,2]); %Add 2nd dimension for pagemtimes

%E is in V/m. Z is in Ohm. B is in T. B/mu = H in A/m
%pagemtimes is a very clever way to do the multiplication that avoids a for
%loop. 10000 times faster.

E = (1/mu)*(pagemtimes(real(Z),real(B))+1i*pagemtimes(imag(Z),imag(B)));
E = permute(E,[1,3,2]);


% parfor i = 1:length(Bx) %Loop through frequencies
%     
%     %E is in V/m. Z is in Ohm. B is in T. B/mu = H in A/m
%     E(:,i) = (1/mu)*(real(Z(:,:,i))*real(B(:,i))+1i*imag(Z(:,:,i))*imag(B(:,i)));
%     
% end
    

Ex = E(1,:);
Ey = E(2,:);
