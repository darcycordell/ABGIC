function [BXpred,BYpred,BZpred,JXext,JYext]=secs(Blat,Blon,BX,BY,BZ,SECSlat,SECSlon,qlat,qlon)
% This code is a slightly modified version of a script from:
% Vanhamaki, H., Juusola, L. (2020). Introduction to spherical
% elementary current systems, in: Dunlop, M., Luhr, H. (eds),
% Ionospheric Multi-Spacecraft Analysis Tools. 
% Doi: https://doi.org/10.1007/978-3-030-26732-2_2
%
% Modified by Darcy Cordell (DC), 2023.
%
%
% Usage: [BXpred,BYpred,BZpred]=secs(Blat,Blon,BX,BY,BZ,SECSlat,SECSlon,qlat,qlon)
%
%   This function takes input magnetic time domain data at N magnetic
%   observatories for M time samples and solves the SECS system for any
%   other set of P data points. The SECS system is then used to solve for
%   the ground magnetic fields at Q query points.
% 
%
% INPUTS:
%       (Blat, Blon): [1 x N] vectors of magnetometer observatory locations
%       BX: [M x N] matrix of X (north) magnetic data in nT
%       BY: [M x N] matrix of Y (east) magnetic data in nT
%       BZ: [M x N] matrix of Z (down) magnetic data in nT
%       (SECSlat, SECSlon): latitude and longitude points where the SECS system 
%           is solved for. These variables can be meshgrid matrices or
%           vectors.
%       (qlat, qlon): latitude and longitude points where the ground magnetic
%           field will be solved for. These variables can be meshgrid
%           matrices or vectors.
%
%       Note: If we are doing a spatial interpolation (as is most common)
%       using a meshgrid of query points (qlat,qlon), then it is most
%       stable to design SECSlat and SECSlon such that they are larger than
%       the area encompassed by the (qlat,qlon) grid. It is also preferable
%       to stagger the (SECSlat,SECSlon) grid from the (qlat,qlon) grid.
%
%       Example: 
%           minlat = 45; maxlat = 60; %the lat boundaries of your desired B grid
%           minlon = -120; maxlon = -105; % the lon boundaries of your desired B grid
%           dlat = 0.5; dlon = 0.5; %grid spacing
%           [qLAT,qLON] = meshgrid(minlat:dlat:maxlat,minlon:dlon:maxlon);
%
%           %Ensure the SECS grid is larger by 2*dlat and 2*dlon and offset
%           % by 1/2 of dlat and lon to stagger the grid:
%           [SECSlat,SECSlon] = meshgrid(minlon-dlat*2-dlat/2:dlat:maxlat+dlat*2+dlat/2,minlon-dlon*2-dlon/2:dlon:maxlon+dlon*2+dlon/2);
%
% OUTPUTS:
%   (BXpred,BYpred,BZpred): [M x Q] matrix of interpolated (i.e. predicted)
%          magnetic field data for N, E and down, respectively, where
%          Q = length(qlat(:))
%
% The original version only output the SECS currents (J), but did not solve
% for the B-fields. So DC modified to include this as output to the
% function. The original version was also written as a stand-alone example
% script, but DC modified it so that it is a more general function
% that can be used with any dataset.
%
%------------------------------ORIGINAL CODE PRE-AMBLE---------------------
%
% This code is free software.
% It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
%
%
%This is an example of using the SECS method to determine the ionospheric
%equivalent current from ground magnetic measurements. 
%  Input: ground magnetometer recordings
%  Output: ionospheric equivalent currents
%
%Based on code written by Ari Viljanen, Antti Pulkkinen, Liisa Juusola and
%Heikki Vanhamaki. Tested with MatLab R2014a. Partly compatible with Octave.
%
%
%
%This code is included as supplementary material to Chapter 2 in the book "Multi-Satellite
%Data Analyses", edited by Malcolm Dunlop and Hermann Luhr. Feel free to use and modify the
%code, but if you use it in any publications, do cite the above mentioned book and the
%following article where the method was originally introduced:
%
%Amm, O. and A. Viljanen, 1999. Ionospheric disturbance magnetic field continuation from
%the ground to the ionosphere using spherical elementary current systems. Earth, Planets
%and Space, 51, 431-440, https://dx.doi.org/10.1186/BF03352247
%
%
%
%  Other references and recommended reading:
%
%Amm, O., 1997. Ionospheric Elementary Current Systems in Spherical Coordinates and Their
%Application. J. Geomag. Geoelectr., 49, 947-955.
%
%Pulkkinen, A., O. Amm, A. Viljanen and BEAR Working Group, 2003a. Ionospheric equivalent
%current distributions determined with the method of spherical elementary current systems.
%J. Geophys. Res., 108, https://dx.doi.org/10.1029/2001JA005085.
%
%Pulkkinen, A., O. Amm, A. Viljanen and BEAR Working Group, 2003b. Separation of the
%geomagnetic variation field on the ground into external and internal parts using the
%spherical elementary current system method, Earth Planets Space, 55, 117–129.
%
%Juusola L., K. Kauristie, H. Vanhamaki, A. Aikio, and M. van de Kamp, 2016: 
%Comparison of Auroral Ionospheric and Field-Aligned Currents Derived From Swarm and
%Ground Magnetic Field Measurements. JGR, https://dx.doi.org/10.1002/2016JA022961
%
%Vanhamaki, H., O. Amm and A. Viljanen, 2003:
%One-dimensional upward continuation of the ground magnetic field disturbance using
%spherical elementary current systems. Earth, Planets and Space, 55, 613-625.
%
%Weygand J., Amm, O., Viljanen, A. et al., 2011. Application and validation of the spherical
%elementary currents systems technique for deriving ionospheric equivalent currents with the
%North American and Greenland ground magnetometer arrays. J. Geophys. Res., 116, A03305,
%doi:10.1029/2010JA016177.
%
%-------------------------END OF ORIGINAL CODE PRE-AMBLE-------------------




%Number of timesteps and stations
[Nt,Nstat] = size(BX);


RE=6371.2;    %Earth's radius [km]
Rext=RE+110;  %radius of the external current layer = ionosphere
epsSVD=0.05;      %epsilon for SVD inversion, see Weygand et al. (2011) for selecting the value 
separaatio=true;  %whether to make the separation into internal and external part or not, see Pulkkinen et al. (2003b)


if separaatio
  %Could and maybe should use different grids for the internal/external SECS and Jeq,
  %but for simplicity use the same. Also we should probably use different epsilon if
  %separation is made, but we also neglect that.
  Rint=RE-20;   %radius of the internal current layer
else
 Rint=[];
end

Neq = length(qlat(:));       % # of Jeq points
Npole = length(SECSlat(:));  % # of SECS poles


%-----------------------------------------
%Determine SECS amplitudes and calculate the equivalent current

%Calculate geometric matrices that give the magnetic field of individual SECS at the magnetometer stations.
%external
[MrE,MtE,MpE] = sub_SECS_2D_DivFree_magnetic((90-Blat)/180*pi,Blon/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,RE,Rext);
if separaatio
  %internal
  [MrI,MtI,MpI] = sub_SECS_2D_DivFree_magnetic((90-Blat)/180*pi,Blon/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,RE,Rint);
else
  MrE=[];  MrI=[];  MtI=[];  MpI=[];
end

%Collect magnetic field components from all stations into a matrix, where each
%column is one timestep and all the components are included in the same column vector	
if separaatio
  %Use all 3 components of B
  Bvecs = NaN(3*Nstat,Nt);
  Bvecs(1:3:3*Nstat,:) = -BZ';  %r=-down
  Bvecs(2:3:3*Nstat,:) = -BX';  %theta=-north
  Bvecs(3:3:3*Nstat,:) =  BY';  %phi=east
else
  %Use horizontal part of B
  Bvecs = NaN(2*Nstat,Nt);
  Bvecs(1:2:2*Nstat,:) = -BX';
  Bvecs(2:2:2*Nstat,:) =  BY';
end

%Make a system matrix that gives the magnetic field components at the stations
%from the SECS scaling factors, Bvecs=matB*Idf, and invert it for Idf
if separaatio
  matB = NaN(3*Nstat,2*Npole);
  matB(1:3:3*Nstat,1:Npole) = MrE;
  matB(2:3:3*Nstat,1:Npole) = MtE;
  matB(3:3:3*Nstat,1:Npole) = MpE;
  matB(1:3:3*Nstat,1+Npole:end) = MrI;
  matB(2:3:3*Nstat,1+Npole:end) = MtI;
  matB(3:3:3*Nstat,1+Npole:end) = MpI;
else
  matB = NaN(2*Nstat,Npole);
  matB(1:2:2*Nstat,:) = MtE;
  matB(2:2:2*Nstat,:) = MpE;
end
%truncated SVD with predefined epsilon
Idf = sub_inv_SVD(matB,epsSVD) * Bvecs;

%Calculate Jeq at the output grid and save also the SECS scaling factors
%external part
[MJtheta,MJphi] = sub_SECS_2D_DivFree_vector((90-qlat(:))/180*pi,qlon(:)/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,Rext,0);
Iext = Idf(1:Npole,:);
JXext = -(MJtheta*Iext)';  %northward component
JYext = (MJphi*Iext)';     %eastward component
if separaatio
  %internal part
  [MJtheta,MJphi] = sub_SECS_2D_DivFree_vector((90-qlat(:))/180*pi,qlon(:)/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,Rint,0);
  Iint = Idf(1+Npole:end,:);
  JXint = -(MJtheta*Iint)';
  JYint = (MJphi*Iint)';
else
  Iint=[];  JXint=[];  JYint=[];
end

%Save the data to a matlab binary file
%Jfile='SECS_theory_example.mat';
%fprintf('\n  Saving data of equivalent currents in %s\n\n',Jfile);
%save(Jfile, 'Blat','Blon','Bnames','BX','BY','BZ','RE','separaatio','Rext','Rint','SECSlat','SECSlon','qlat','qlon','Iext','Iint','JXext','JYext','JXint','JYint');

%%
[PBrE,PBtE,PBpE] = sub_SECS_2D_DivFree_magnetic((90-qlat)/180*pi,qlon/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,RE,Rext);
if separaatio
  %internal
  [PBrI,PBtI,PBpI] = sub_SECS_2D_DivFree_magnetic((90-qlat)/180*pi,qlon/180*pi,(90-SECSlat(:))/180*pi,SECSlon(:)/180*pi,RE,Rint);
else
  PBrE=[];  PBrI=[];  PBtI=[];  PBpI=[];
end

%Make a system matrix that gives the magnetic field components at the stations
%from the SECS scaling factors, Bvecs=matB*Idf, and invert it for Idf
if separaatio
  matBP = NaN(3*Neq,2*Npole);
  matBP(1:3:3*Neq,1:Npole) = PBrE;
  matBP(2:3:3*Neq,1:Npole) = PBtE;
  matBP(3:3:3*Neq,1:Npole) = PBpE;
  matBP(1:3:3*Neq,1+Npole:end) = PBrI;
  matBP(2:3:3*Neq,1+Npole:end) = PBtI;
  matBP(3:3:3*Neq,1+Npole:end) = PBpI;
else
  matBP = NaN(2*Neq,Npole);
  matBP(1:2:2*Neq,:) = PBtE;
  matBP(2:2:2*Neq,:) = PBpE;
end

Bpred = mtimes(matBP,Idf);

BXpred = -Bpred(2:3:3*Neq,:);
BYpred = Bpred(3:3:3*Neq,:);
BZpred = -Bpred(1:3:3*Neq,:);


%############################################## END OF FUNCTION SECS_theory_example


function invM=sub_inv_SVD(M,epsSVD)
%
%Calculate the inverse of a given matrix M using
%singular value decomposition (SVD).
%
%       M: Original matrix
%  epsSVD: Regularization parameter, scalar, 0 <= epsSVD <= 1
%    invM: Inverse matrix
%
% Functions called:
% - sub_where


l1=length(M(:,1));
l2=length(M(1,:));

fprintf('\n  Calculating SVD of a [%d,%d] matrix ...',l1,l2);
[U,S,V]=svd(M,'econ');
fprintf(' done\n');

lkms=min(l1,l2);
%Vector s is ordered so that s(n)>=s(n+1) and all s(n)>=0.
s=diag(S(1:lkms,1:lkms));
U=U(:,1:lkms);
V=V(:,1:lkms);

%Calculate the inverse matrix
slim=epsSVD*s(1);
lkm1=sum(s<=slim);
fprintf('  epsilon = %f, singular values range from %f to %f \n', epsSVD, s(1), s(lkms));
fprintf('  --> %d values smaller than %f deleted (of %d values)\n\n', lkm1, slim, lkms);
invM = V * diag(sub_where(s<=slim,0,1./s)) * U';


%############################################## END OF FUNCTION sub_inv_SVD


function [y]=sub_where(c,a,b)
%
%TeLa-style where-function for matlab.
%
%[y] = where(c,a,b)
% where(c,a,b) returns a where c is nonzero and b where c is zero.
%   Usually the arguments are arrays, and if they are, they must be
%   of similar dimensions.

%convert scalar a or b to array
a=a+0*c;
b=b+0*c;

y=a;
ind=find(c==0);
y(ind)=b(ind);


%############################################## END OF FUNCTION sub_where


function [matBradial,matBtheta,matBphi]=sub_SECS_2D_DivFree_magnetic(thetaB,phiB,thetaSECS,phiSECS,Rb,Rsecs)
%
%Matlab function for calculating matrices matB? that relate the (ground) magnetic field to the
%SECS representation so that
%
%  B? = mat? * Idf, where
%   ? = {r,theta,phi}             vector component
%  B? = [B?(1) B?(2) ... ]'       measurements
%  Idf = [Idf(1') Idf(2') ... ]'  scaling factors
%
%Assume units: [B]=nT, [Idf]=A and [length]=km
%
%
%       thetaB,phiB : Co-latitude and longitude of the magnetometer stations, [radian], Nb-dimensional vectors
% thetaSECS,phiSECS : Co-latitude and longitude of the DF SECS, [radian], Nsecs-dimensional vectors
%                Rb : Radial distance of the magnetometer stations [km], scalar or Nb-dimensional vector
%          Rb,Rsecs : Radii of the spheres where magnetometers and DF SECS are located, [km], scalars
%
% matBradial,matBtheta,matBphi : (Nb,Nsecs)-matrices that relate SECS scaling factors to the magnetic field
%
%
%Heikki V, March 2010, matlab 7.8.0.347 (R2009a)


%number of points where B is calculated and scaling factors are given
Nb=length(thetaB(:));
Nsecs=length(thetaSECS(:));
matBradial=NaN(Nb,Nsecs);
matBtheta=NaN(Nb,Nsecs);
matBphi=NaN(Nb,Nsecs);

if length(Rsecs)~=1
  fprintf('\n  Sorry, input parameter Rsecs must be scalar \n\n');
  return;
end

%If Rb is scalar, use same radius for all points.
Rb=Rb+0*thetaB;

%Ratio of the radii, smaller/larger
suhde=min(Rb,Rsecs) ./ max(Rb,Rsecs);

%There is a common factor mu0/(4*pi)=1e-7. Also 1/Rb is a common factor
%If scaling factors are in [A], radii in [km] and magnetic field in [nT]  --> extra factor of 1e6
kerroin=0.1 ./ Rb;

%loop over B positions
for n=1:Nb
  %cos and square of sin of co-latitude in the SECS-centered system
  %See Eq. (A5) and Fig. 14 of Vanhamaki et al.(2003)
  CosThetaPrime=cos(thetaB(n))*cos(thetaSECS) + sin(thetaB(n))*sin(thetaSECS).*cos(phiSECS-phiB(n));
  Sin2ThetaPrime=(1-CosThetaPrime.^2);

  %sin and cos of angle C, divided by sin(theta').
  %See Eqs. (A2)-(A5) and Fig. 14 of Vanhamaki et al.(2003)
  ind=find((Sin2ThetaPrime>1e-10));
  SinC=zeros(size(CosThetaPrime));
  CosC=zeros(size(CosThetaPrime));
  SinC(ind) = sin(thetaSECS(ind)).*sin(phiSECS(ind)-phiB(n)) ./ Sin2ThetaPrime(ind);
  CosC(ind) = (cos(thetaSECS(ind)) - cos(thetaB(n))*CosThetaPrime(ind)) ./ (sin(thetaB(n))*Sin2ThetaPrime(ind));

  %auxiliary variable
  juuri=sqrt(1 - 2*suhde(n)*CosThetaPrime + suhde(n)^2);

  if Rb(n) < Rsecs
    apuVertical = 1;
    %See Eq. (10) of Amm and Viljanen (1999).
    apuHorizontal = -kerroin(n) * ((suhde(n)-CosThetaPrime)./juuri + CosThetaPrime);
  elseif Rb(n) > Rsecs
    apuVertical = suhde(n);
    %See Eq. (A8) of Amm and Viljanen (1999).
    apuHorizontal = -kerroin(n) * ((1-suhde(n)*CosThetaPrime)./juuri - 1);
  else
    %Rb(n) == Rsecs
    apuVertical = 1;
    %Actually horizontal field is not well defined, but this is the average.
    %See Eqs. (10) and (A8) of Amm and Viljanen (1999).
    apuHorizontal = -kerroin(n) * (juuri+CosThetaPrime-1)/2;
  end

  matBradial(n,:) = apuVertical * kerroin(n) * (1./juuri(:) - 1);  %See Eqs. (9) and (A7) of Amm and Viljanen (1999).
  matBtheta(n,:) = apuHorizontal(:) .* CosC(:);
  matBphi(n,:) = -apuHorizontal(:) .* SinC(:);
end


%############################################## END OF FUNCTION sub_SECS_2D_DivFree_magnetic


function [matVtheta,matVphi]=sub_SECS_2D_DivFree_vector(thetaV,phiV,thetaSECS,phiSECS,radius,LimitAngle)
%
%Matlab function for calculating matrices matVtheta and matVphi which give the theta- and
%phi-components of a vector field from the scaling factors of div-free spherical elementary
%current systems (DF SECS),
%  Vtheta = matVtheta * Idf, where
%  Vtheta = [Vtheta(1) Vtheta(2) ...]'  vector of theta-components
%  Idf = [Idf(1') Idf(2') ...]' vector of scaling factors
%
%
%  INPUT
%        thetaV,phiV : Co-latitude and longitude of points where the vector field is to be calculated,
%                       [radian], Nv-dimensional vectors
%  thetaSECS,phiSECS : Co-latitude and longitude of the DF SECS, [radian], Nsecs-dimensional vectors
%             radius : Radius of the sphere where the calculation takes place, [km], scalar
%         LimitAngle : Half-width of the uniformly distributed SECS, [radian], scalar or Nsecs-dimensional vector
%
%  OUTPUT
%  matVtheta,matVphi : (Nv,Nsecs)-matrices that relate SECS scaling factors to the vector field
%
%
%  NOTE: Each individual SECS is assumed to be uniformly distributed inside a spherical cap with
%        half angle 'LimitAngle'. This removes the singularity at the pole of the SECS. Outside the
%        cap this kin of SECS has exactly the same field as the traditional singular SECS defined
%        by Amm (1998).
%
%
%Heikki V, March 2010, matlab 7.8.0.347 (R2009a)


%theta- and phi-directions are badly defined at the poles of the spherical coordinate
%system, but no check is done here...


%number of points where V is calculated and scaling factors are given
Nv=length(thetaV(:));
Nsecs=length(thetaSECS(:));
matVtheta=NaN(Nv,Nsecs);
matVphi=NaN(Nv,Nsecs);

%if LimitAngle is scalar, use it for every SECS
if length(LimitAngle(:))==1
  LimitAngle=LimitAngle+zeros(size(thetaSECS));
end;

%This is a common factor in all components
CommonFactor=1/(4*pi*radius);

%loop over vector field positions
for n=1:Nv
  %cosine of co-latitude in the SECS-centered system
  %See Eq. (A5) and Fig. 14 of Vanhamaki et al.(2003)
  CosThetaPrime=cos(thetaV(n))*cos(thetaSECS) + sin(thetaV(n))*sin(thetaSECS).*cos(phiSECS-phiV(n));

  %sin and cos of angle C, multiplied by sin(theta').
  %See Eqs. (A2)-(A5) and Fig. 14 of Vanhamaki et al.(2003)
  SinC = sin(thetaSECS).*sin(phiSECS-phiV(n));
  CosC = (cos(thetaSECS) - cos(thetaV(n))*CosThetaPrime) / sin(thetaV(n));

  %Find those SECS poles that are far away from the calculation point
  distant=(CosThetaPrime < cos(LimitAngle));

  %vector field proportional to cot(0.5*CosThetaPrime), see Eq. (2) of Vanhamaki et al.(2003)
  dummy = CommonFactor ./ (1-CosThetaPrime(distant));
  matVtheta(n,distant) = dummy .* SinC(distant);
  matVphi(n,distant) = dummy .* CosC(distant);

  %Assume that the curl of a DF SECS is uniformly distributed inside LimitAngle
  % --> field proportional to a*tan(0.5*CosThetaPrime), where a=cot(0.5*LimitAngle)^2
  dummy = CommonFactor * cot(0.5*LimitAngle(~distant)).^2 ./ (1+CosThetaPrime(~distant));
  matVtheta(n,~distant) = dummy .* SinC(~distant);
  matVphi(n,~distant) = dummy .* CosC(~distant);
end;

