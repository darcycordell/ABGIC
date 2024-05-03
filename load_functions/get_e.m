function [ex,ey] = get_e(Z,Zf,Bx,By,b)

disp('*******************CALCULATE E******************************')

fAxis = b(1).fAxis;
f = b(1).f;
nf = length(b(1).fAxis);
nt = b(1).nt;
ns = size(Z,3);

%Compute E-field using 3-D impedances
Ex = nan(nf,ns)+1i*zeros(nf,ns); Ey = nan(size(Ex));
ex = nan(nt,ns); ey = nan(nt,ns);
for i = 1:ns
    Zint3D =  interpolate_Z(Z(:,:,i),Zf,f,fAxis);
    [Ex(:,i),Ey(:,i)] = calc_E(Bx(i,:),By(i,:),Zint3D);

    [ex(:,i),ey(:,i)] = calc_ifft(Ex(:,i),Ey(:,i),b(1).pad);

    disp(['E field for MT Site #: ',num2str(i),' completed'])
end

disp('*******************E: Completed******************************')




