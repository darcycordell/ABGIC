function [a,zn,fwd] = calc_Z_trich(d,zn,in)
%
% Function which sets up the 1-D resistivity models from Trichtchenko et al. (2019)
% and calculates the 1-D impedance for each model
%
% Inputs: 
%       d = MT impedance data structure
%       zn = Zone geographical information
%       in = logical mask telling which MT site belongs in which zone
%
% Outputs:
%       a = data structure identical to "d" but with 1-D impedances at MT
%       site locations
%       zn = Updated zone information with models, depths, etc.
%       fwd = raw forward modelled MT data for each zone

%Zone #7 (from Trichtchenko et al. 2017)
zn(1).thick = [50 1400 11600 16000 9000 62000 150000 160000 110000 150000 230000 100000];
zn(1).rho = [30 10 3000 1750 1200 1400 635 50 20 5.62 1.58 1.12];

%Unnamed Zone
zn(2).thick = zn(1).thick; %Zone #2 (zn(2)) has no resistivity model associated with it but it is              
zn(2).rho = zn(1).rho;     %given the same values as Zone #1 here

%Zone #6
zn(3).thick = [65 2700 9000 17000 10000 61000 150000 160000 110000 150000 230000 100000];
zn(3).rho = [50 15 3000 2300 2100 1400 160 50 20 5.6 1.58 1.12];

%Zone #5
zn(4).thick = [35 3100 16000 10000 10000 61000 150000 160000 110000 150000 230000 100000];
zn(4).rho = [30 20 440 900 330 950 160 40 20 5.6 1.58 0.89];

%Zone #10
zn(5).thick = [175 1900 6000 10000 22000 60000 150000 160000 110000 150000 230000 100000];
zn(5).rho = [30 10 1000 1300 800 500 625 50 20 5.62 1.58 1.12];

%Zone #8
zn(6).thick = [25 3300 11700 15000 9000 61000 150000 160000 110000 150000 230000 100000];
zn(6).rho = [30 10 900 275 360 315 660 50 20 5.62 1.58 1.12];

%Zone #3
zn(7).thick = [40 2900 16000 10000 12000 59000 150000 160000 110000 150000 230000 100000];
zn(7).rho = [30 6 2000 3000 620 150 230 90 8 2.4 0.89 0.47];

%Zone #1
zn(8).thick = [75 2200 13000 12500 17500 55000 150000 160000 110000 150000 230000 100000];
zn(8).rho = [30 30 385 2500 2500 2250 550 40 8 2.4 0.89 0.47];

%Zone #9
zn(9).thick = [25 1900 10100 12000 16000 60000 150000 160000 110000 150000 230000 100000];
zn(9).rho = [15 10 4200 4500 4500 1600 680 50 20 5.62 1.58 1.12];

%Zone #2
zn(10).thick = [50 2400 21600 13000 8000 55000 150000 160000 110000 150000 230000 100000];
zn(10).rho = [30 25 3000 2000 4000 2000 260 70 8 2.4 0.89 0.47];

%Zone #4
zn(11).thick = [50 2000 17000 10000 14000 57000 150000 160000 110000 150000 230000 100000];
zn(11).rho = [50 10 1900 1700 850 500 530 90 8 2.4 0.89 0.47];

%Note that Zones 12, 13, 14, 15, 17, 18, 19, 20, and 24 do not have MT sites in them
%but I've included them here in case I ever decide to do a 3-D model
%comparison of the region

% Zone #1
zn(12).thick = [25 1000 8500 10500 10000 70000 50000 100000 160000 110000 150000 230000 100000];
zn(12).rho = [30 25 610 310 90 200 230 200 40 11 2 1.22 0.77];

% Zone 2N
zn(13).thick = [2 12500 15500 3500 68500 150000 160000 110000 150000 230000 100000];
zn(13).rho = [115 43000 34800 26000 25800 2300 100 20 5.62 1.58 0.89];

% Zone 3B
zn(14).thick = [50 2500 8500 12000 9000 68000 50000 100000 160000 110000 150000 230000 100000];
zn(14).rho = [75 50 600 550 415 210 95 380 80 25 5.62 1.58 0.89];

% Zone 5N
zn(15).thick = [10 11500 16000 11500 61000 150000 160000 110000 150000 230000 100000];
zn(15).rho = [70 43800 24900 46200 1750 1950 120 20 5.62 1.58 0.89];

%Zone 4S
zn(16).thick = [2 14500 9500 13000 63000 50000 100000 160000 110000 150000 230000 100000];
zn(16).rho = [70 740 300 65 130 265 210 50 20 5.62 1.58 0.89];

% Zone 3N
zn(17).thick = [10 13000 20000 3500 63500 150000 160000 110000 150000 230000 100000];
zn(17).rho = [70 27700 30900 21000 15600 3200 100 20 5.62 1.58 0.89];

% Zone 1
zn(18).thick = [25 1000 8500 10500 10000 70000 50000 100000 160000 110000 150000 230000 100000];
zn(18).rho = [30 25 610 310 90 200 230 200 40 11 2 1.22 0.77];

% Zone 2S
zn(19).thick = [2 11000 12000 11500 65500 50000 100000 160000 110000 150000 230000 100000];
zn(19).rho = [115 820 810 260 120 120 210 50 20 5.62 1.58 0.89];

% Zone 4N
zn(20).thick = [15 16500 12000 9000 62500 150000 160000 110000 150000 230000 100000];
zn(20).rho = [70 27800 1400 20 16350 1850 100 20 5.62 1.58 0.89];

%Zone 5S
zn(21).thick = [2 13000 7500 21500 58000 50000 100000 160000 110000 150000 230000 100000];
zn(21).rho = [70 520 500 250 370 750 210 50 20 5.62 1.58 0.89];

% Zone 6
zn(22).thick = [50 3800 9200 6000 20000 61000 150000 160000 110000 150000 230000 100000];
zn(22).rho = [70 100 14000 7700 12000 3300 700 120 20 2.4 1.12 0.48];

% Zone 3S
zn(23).thick = [10 14000 11000 8000 67000 50000 100000 160000 110000 150000 230000 100000];
zn(23).rho = [70 850 400 300 140 80 210 50 20 5.62 1.58 0.89];

% Zone 2V
zn(24).thick = [250 2500 8250 9000 15000 65000 50000 100000 160000 110000 150000 230000 100000];
zn(24).rho = [100 25 75 150 110 115 150 210 50 20 5.62 1.58 0.89];

% % Saskatchewan Zones (for 3-D modelling in the padding)
% 
% % Zone 1A
% zn(25).thick = [3 13500 13500 11000 62000 150000 160000 110000 150000 230000 100000];
% zn(25).rho = [100 16900 4350 4350 500 200 29 8 2.4 1.12 0.48];
% 
% % Zone 1B
% zn(26).thick = [100 1700 14200 12000 9000 63000 150000 160000 110000 150000 230000 100000];
% zn(26).rho = [25 10 2300 2500 2000 1600 180 29 8 2.4 1.12 0.48];
% 
% % Zone 1C
% zn(27).thick = [15 1450 12000 13500 9500 63000 150000 160000 110000 150000 230000 100000];
% zn(27).rho = [100 2200 12800 4350 4350 1600 180 29 8 2.4 1.12 0.48];
% 
% % Zone 1D
% zn(28).thick = [100 1700 14200 13000 13000 58000 150000 160000 110000 150000 230000 100000];
% zn(28).rho = [25 10 2200 1950 900 550 350 90 8 2.4 1.12 0.48];
% 
% % Zone 1E
% zn(29).thick = [50 2400 19600 13000 11000 54000 150000 160000 110000 150000 230000 100000];
% zn(29).rho = [25 10 1700 3200 4000 4000 425 70 8 2.4 1.12 0.48];
% 
% % Zone 2A
% zn(30).thick = [3 13500 14500 10000 62000 150000 160000 110000 150000 230000 100000];
% zn(30).rho = [115 2500 2300 900 600 200 29 8 2.4 1.12 0.48];
% 
% % Zone 2B
% zn(31).thick = [15 400 13100 14500 10000 62000 150000 160000 110000 150000 230000 100000];
% zn(31).rho = [115 2200 4600 2300 900 600 200 29 8 2.4 1.12 0.48];
% 
% %Zone 3A
% zn(32).thick = [5 18000 18000 64000 150000 160000 110000 150000 230000 100000];
% zn(32).rho = [100 37500 11500 23600 24800 500 8 2.4 1.12 0.48];
% 
% % Zone 3B
% zn(33).thick = [100 200 17700 18000 64000 150000 160000 110000 150000 230000 100000];
% zn(33).rho = [100 10 37500 11500 23600 24800 500 8 2.4 1.12 0.48];
% 
% % Zone 3C
% zn(33).thick = [25 750 17250 18000 64000 150000 160000 110000 150000 230000 100000];
% zn(33).rho = [100 2200 37500 11500 23600 24800 500 8 2.4 1.12 0.48];
% 
% % Zone 4A
% zn(34).thick = [2 16000 12000 12000 60000 150000 160000 110000 150000 230000 100000];
% zn(34).rho = [100 9200 12500 5800 450 160 29 8 2.4 1.12 0.48];
% 
% % Zone 4B
% zn(34).thick = [100 1400 11000 26500 7000 54000 150000 160000 110000 150000 230000 100000];
% zn(34).rho = [25 5 350 375 300 500 300 29 8 2.4 1.12 0.48];
% 
% % Zone 4F
% zn(34).thick = [1 14500 9000 20500 56000 150000 160000 110000 150000 230000 100000];
% zn(34).rho = [70 1700 1500 1700 500 500 29 8 2.4 1.12 0.48];
% 
% % Zone SBZb
% zn(34).thick = [75 1300 11200 22500 10000 55000 150000 160000 110000 150000 230000 100000];
% zn(34).rho = [25 5 500 100 100 80 55 29 8 2.4 1.12 0.48];

nb = length(zn);
for i = 1:nb
    zn(i).depth = [0 cumsum(zn(i).thick)];
    zn(i).depth = zn(i).depth(1:end-1);
end


% Compute 1-D MT curves for each 1-D Zone Model and assign this data to a 1-D data structure
for i = 1:nb
    [fwd(i)]=calc_fwd_1d(zn(i).depth,zn(i).rho,d.f',0);
end

a = d;

a.Z(:,[1 4],:) = 0 +1i*0;

for is = 1:a.ns
    id = find(in(is,:)==1);
    a.Z(:,2,is) = fwd(id).Z;
    a.Z(:,3,is) = -fwd(id).Z;
end

[a.rho,a.pha,a.rhoerr,a.phaerr] = calc_rho_pha(a.Z,a.Zerr,a.T);

ind = isnan(d.Z);
a.Z(ind) = NaN+1i*NaN;
a.rho(ind) = NaN;
a.pha(ind) = NaN;
a.Zerr(ind) = NaN+1i*NaN;
a.rhoerr(ind) = NaN;
a.phaerr(ind) = NaN;