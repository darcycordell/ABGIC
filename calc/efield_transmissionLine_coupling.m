function [avg_strike,avg_r,std_strike,emag,simple_LV] = efield_transmissionLine_coupling(exnn3d,eynn3d,dx,dy,d)

line_angle = (180/pi)*mod(atan2(dy,dx),2*pi);



e_angle = (180/pi)*mod(atan2(eynn3d,exnn3d),2*pi);
emag_all = sqrt(eynn3d.^2+exnn3d.^2);
emag = mean(emag_all,1);
coupling_angle = line_angle'-e_angle;

strike = coupling_angle';
J = 1;

strike=strike*J;
strike=strike*(pi/180); %convert to radians
n=length(strike(:,1));

%Take the cosine and sine of the strike angles (the x and y projections)
cosine=cos(strike);
sine=sin(strike);

% average the cos and sin projections
avg_c=sum(cosine,2)/n;
avg_s=sum(sine,2)/n;

%The resultant is 0 < r < 1 and is a measure of scatter. As you add the
%vectors, if they all pointed the same way, you would have r=1. If you had
%all the vectors pointing different directions, then r=0.
avg_r=sqrt(avg_c.^2+avg_s.^2);

%The average radial strike
avg_strike=atan2(avg_s,avg_c)*(180/pi)/J;

%The average radial standard deviation in degrees.
std_strike=(1/J)*sqrt((-2*log(avg_r)))*(180/pi);


simple_LV = d*emag.*cosd(avg_strike');