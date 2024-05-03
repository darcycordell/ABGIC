function d = get_Z(option_flag)

get_files;

curdir = pwd;

if option_flag == 1 %Real 3-D Impedance data loaded from mat file on path
    
    cd(zFilePath)
    load(zFile);
    cd(curdir)
    disp(['Loaded MT Impedance File: ',zFile])


elseif option_flag == 2 %Uniform halfspace impedance everywhere

    [L,~,~] = get_network;

    rx = []; ry = [];
    for sidx = 1:length(L)
                
        rx = [rx L(sidx).Loc(1:6:end,1)'];
        ry = [ry L(sidx).Loc(1:6:end,2)'];
       
    end

    d.loc(:,1) = rx'; d.loc(:,2) = ry';

    d.loc = unique(d.loc,'rows');

    d.ns = length(d.loc(:,1));

    for i = 1:d.ns
        d.site{i,1} = ['Site_',sprintf('%03.f',i)];
    end

    
    d.nf = 40;
    d.f = 10.^linspace(1,-4.5,d.nf)';

    d.Z = zeros(d.nf,4,d.ns);

    mu0 = 4*pi*10^-7;
    sigma = 0.01;

    for i = 1:d.ns
        d.Z(:,2,i) = ((1+1i)/sqrt(2))*sqrt((2*pi*d.f'*mu0)./sigma);
        d.Z(:,3,i) = -((1+1i)/sqrt(2))*sqrt((2*pi*d.f'*mu0)./sigma);
    end

    d.Zerr = nan(size(d.Z));
    [d.rho,d.pha] = calc_rho_pha(d.Z,d.Zerr,1./d.f);

    disp(['Using Impedance for halfspace with value: ',num2str(1./sigma),' Ωm'])

elseif option_flag == 3 %Trichtchenko 1-D Piecewise Zones

    [zn, Clat, Clon, txt, nb] = load_trich_zones;

    cd(zFilePath)
    [d,in] = load_assign_impedance(zn,zFile);

    d.Z(abs(real(d.Z(:)))>10^5)=NaN;

    [d] = calc_Z_trich(d,zn,in);

end



