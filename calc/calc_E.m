function [Ex,Ey,ZxxBx,ZxyBy,ZyxBx,ZyyBy] = calc_E(Bx,By,Z,f)
% Function to calculate the E-field vector given B-field vector and 2x2
% impedance in frequency domain
%
% Inputs: Bx and By are north (x) and east (y) magnetic fields in nT at a
%           specific location in the frequency domain
%           Size is (1 x nf)
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

E = (1/mu)*pagemtimes(Z,B);
E = permute(E,[1,3,2]);


% parfor i = 1:length(Bx) %Loop through frequencies
%     
%     %E is in V/m. Z is in Ohm. B is in T. B/mu = H in A/m
%     E(:,i) = (1/mu)*(real(Z(:,:,i))*real(B(:,i))+1i*imag(Z(:,:,i))*imag(B(:,i)));
%     
% end

%ZxxBx = (1/mu)*squeeze(Z(1,1,:).*B(1,:,:));
%ZxyBy = (1/mu)*squeeze(Z(1,2,:).*B(2,:,:));
%ZyxBx = (1/mu)*squeeze(Z(2,1,:).*B(1,:,:));
%ZyyBy = (1/mu)*squeeze(Z(2,2,:).*B(2,:,:));
    

Ex = E(1,:);
Ey = E(2,:);

%%
% loglog(1./flip(f(f>0)),squeeze(abs(Z(1,2,f>0)))); hold on;
% loglog(1./flip(f(f>0)),abs(By(f>0))*10^-9)
% loglog(1./flip(f(f>0)),abs(Ex(f>0)))
