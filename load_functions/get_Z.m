function d = get_Z(option_flag)

get_files;

curdir = pwd;

if option_flag == 1 %Real 3-D Impedance data loaded from mat file on path
    
    cd(zFilePath)
    load(zFile);
    cd(curdir)
    disp(['Loaded MT Impedance File: ',zFile])


elseif option_flag == 2 %Uniform halfspace impedance everywhere

    cd(zFilePath)
    load(zFile);
    cd(curdir)

    mu0 = 4*pi*10^-7;
    sigma = 0.01;

    d.Z(:,[1 4],:) = 0+0*1i;
    for i = 1:d.ns
        d.Z(:,2,i) = ((1+1i)/sqrt(2))*sqrt((2*pi*d.f'*mu0)./sigma);
        d.Z(:,3,i) = -((1+1i)/sqrt(2))*sqrt((2*pi*d.f'*mu0)./sigma);
    end

    d.Zerr = nan(size(d.Z));
    [d.rho,d.pha] = calc_rho_pha(d.Z,d.Zerr,1./d.f);

    disp(['Using Impedance for halfspace with value: ',num2str(1./sigma),' Ωm'])

elseif option_flag == 3 %Trichtchenko 1-D Piecewise Zones

    [zn] = load_trich_zones;

    cd(zFilePath)
    [d,in] = load_assign_impedance(zn,zFile);

    d.Z(abs(real(d.Z(:)))>10^5)=NaN;

    [d,zn] = calc_Z_trich(d,zn,in);
    d.in = in;
    d.zn = zn;

elseif option_flag == 4 %Data space method

    cd('/Users/darcycordell/Documents/GitHub/ABGIC/000_DATA/05_MT_IMPEDANCE')
    load('AB_BC_MT_DATA_526_sites_230321.mat');
    cd(curdir)
    disp(['Loaded MT Impedance File: ',zFile])

elseif option_flag == 5 %1-D Model

    cd(zFilePath)
    load(zFile);
    cd(curdir)

    min_freq = 0.0001;
    max_freq = 10;
    num_freq = 80;
    
    d.f = 10.^linspace(log10(min_freq),log10(max_freq),num_freq).';
    d.T = 1./d.f;
    d.nf = length(d.f);
    d.Z = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns);
    
    model_depth = [0 3300];
    model_res = [10 10000];
    %model_res = [10000 10];
    
    
    [fwd]=calc_fwd_1d(model_depth,model_res,d.f',0);
    
    Z = fwd.Z.';

    d.Z(:,[1 4],:) = 0+0*1i;
    for i = 1:d.ns
        d.Z(:,2,i) = Z;
        d.Z(:,3,i) = -Z;
    end

    d.Zerr = nan(size(d.Z));
    [d.rho,d.pha] = calc_rho_pha(d.Z,d.Zerr,1./d.f);

elseif option_flag == 6 %2-D MARE2DEM MODEL

    cd(zFilePath)
    load(zFile);
    cd(curdir)

    cd('/Users/darcycordell/Library/CloudStorage/GoogleDrive-dcordell@ualberta.ca/My Drive/Work/Projects/GIC/Horton_Synthetics/MARE2DEM/run03')
    
    %Model 2 (100 Ωm coast)
    st = m2d_readEMData2DFile('fwd_ocean_500km_5f1000.dat.0.resp');
    sta = m2d_readEMData2DFile('fwd_ocean_500km_5fMT1000.dat.0.resp');

    cd(curdir);
    
    f = st.stMT.frequencies';
    y = st.stMT.receivers(:,2)/1000;
    nf = length(f);
    ns = length(st.stMT.receiverName);

    d.ns = 1;
    d.f = f;
    d.T = 1./d.f;
    d.nf = length(d.f);
    d.Z = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns);
    
    %Load rho/pha
    rhoxy = reshape(10.^sta.DATA(sta.DATA(:,1)==123,7),ns,nf);
    rhoyx = reshape(10.^sta.DATA(sta.DATA(:,1)==125,7),ns,nf);
    phaxy = reshape(sta.DATA(sta.DATA(:,1)==104,7),ns,nf);
    phayx = reshape(sta.DATA(sta.DATA(:,1)==106,7),ns,nf);

    mu0 = 4*pi*10^-7;

    indy = nearestpoint(-5,y);
    d.Z(:,[1 4],:) = 0+0*1i;
    for i = 1:d.ns
        d.Z(:,2,i) = sqrt(mu0*2*pi*f.*rhoxy(indy,:)).*exp(1i*phaxy(indy,:)*pi/180);
        d.Z(:,3,i) = sqrt(mu0*2*pi*f.*rhoyx(indy,:)).*exp(1i*(phayx(indy,:)-180)*pi/180);
    end

    d.Zerr = nan(size(d.Z));
    [d.rho,d.pha] = calc_rho_pha(d.Z,d.Zerr,1./d.f);

end

cd(curdir);


