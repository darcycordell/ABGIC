clearvars

get_files;

cd(mainPath);
fullTic = tic;

%% GEOPHYSICS STEP

if ~uniform

[d] = get_Z(1); %Option flag 1 = MT impedance (model space or data space depending on get_files)
                %Option flag 2 = 1000 Ωm halfspace
                %Option flag 3 = 1-D Trichtchenko
                %Option flag 4 = Data space method hard-coded
%[d1] = get_Z(1);
%d1.Z = d1.Z*1.1;
%[d] = get_Z(4); 


[GIC_Hall,tgic] = get_hall_data;
tgic_orig = tgic;


%Compute E field at each impedance query point
latq = d.loc(:,1);
lonq = d.loc(:,2);


b = get_b('mat');

disp('Interpolating b(t) Spatially...')

[bx,by] = interpolate_b(b,latq,lonq);
disp('...done interpolation')
disp('...converting interpolated b to frequency domain ...')

ttic = tic;
nf = length(b(1).fAxis);
Bx = zeros(length(latq),nf); Bx = Bx+1i*Bx;
By = Bx;
for i = 1:size(bx,1)
    [Bx(i,:),By(i,:)] = calc_fft(bx(i,:)',by(i,:)',b(1).fs,b(1).pad); 
end


disp(['...completed B(omega) in ',num2str(toc(ttic)),' seconds'])

[ex,ey] = get_e(d.Z,d.f,Bx,By,b);

%[d1] = get_Z(3);
% [ex1,ey1] = get_e(d1.Z,d.f,Bx,By,b);
% mage1 = sqrt(ex1.^2+ey1.^2);
% meanE1 = mean(mage1,2,'omitnan');

mage = sqrt(ex.^2+ey.^2);
meanE = mean(mage,2,'omitnan');
[maxmeanE,indmax] = max(meanE);
strike = (180/pi)*atan2(ey,ex);
maxstrike = strike(indmax);


c=cos(strike*pi/180);
s=sin(strike*pi/180);
avg_c=sum(c,2)/size(strike,2);
avg_s=sum(s,2)/size(strike,2);
avg_r=sqrt(avg_c.^2+avg_s.^2);

%The average radial strike
avg_strike=atan2(avg_s,avg_c)*(180/pi);


xmax = repmat(meanE.*cosd(avg_strike),1,d.ns);
ymax = repmat(meanE.*sind(avg_strike),1,d.ns);


[maxe, maxeIndices] = max(mage);
[maxeCounts, maxeGroups] = groupcounts(maxeIndices');
tind = maxeGroups(find(maxeCounts==max(maxeCounts)));


clear Bx By

%Used for Cordell 2024 SW paper
%t1 = datetime('2023-04-24 01:00:00');
%t2 = datetime('2023-04-24 09:00:00');

%t1 = datetime('2023-04-24 2:30:40');
%t2 = datetime('2023-04-24 2:30:50');

%t1 = datetime('2012-03-09 00:00:00');
%t2 = datetime('2012-03-09 20:00:00');

%t1 = datetime('2024-05-11 06:00:00');
%t2 = datetime('2024-05-11 12:00:00');
%t2 = datetime('2024-05-11 07:00:00');

%t1 = datetime('2024-05-11 8:00:00');
%t2 = datetime('2024-05-11 17:00:00');

%t1 = datetime('2021-11-04 06:00:00');
%t2 = datetime('2021-11-04 13:00:00');

%t1 = datetime('2023-03-23 04:00:00');
%t2 = datetime('2023-03-23 18:00:00');

%t1 = datetime('2023-03-24 01:00:00');
%t2 = datetime('2023-03-24 08:00:00');

%t1 = datetime('2023-04-24 7:00:00');
%t2 = datetime('2023-04-24 7:40:00');

%t1 = datetime('2024-9-17 6:00:00');
%t2 = datetime('2024-9-17 14:00:00');

t1 = datetime('2024-10-10 14:45:00');
t2 = datetime('2024-10-10 15:45:00');


tind = find(isbetween(b(1).times,t1,t2));

%tind = find(b(1).times=='2023-04-24 07:19:50');

%tind = 1:length(b(1).times);
%t1 = b(1).times(1); t2 = b(1).times(end);

GIC_Times = b(1).times(tind);
GIC_Times_orig = GIC_Times;

end
%% ENGINEERING STEP
%clearvars
%Load network


%"L" is a nLinesx1 structure with following variables:
%   line name, line voltage, line start point
%   substation, line end point substation, raw line coordinates, and from
%   bus and to bus indices
%"S" is a nSubx1 structure containing substation names, locations, owner and
%earthing resistance
% "T" is a nTransx1 structure containing W1 (HV or Series) winding
% resistance, W2 (LV or Common) winding resistance, HV type, LV type, and
% from bus and to bus indices
[L,S,T] = get_network;
%L((find(contains({L(:).Name},'1057')))) = [];

%L(100).Resistance = L(100).Resistance*0.3;
%L(106).Resistance = L(106).Resistance*0.3;
%L(108).Resistance = L(108).Resistance*0.3;
%L(114).Resistance = L(114).Resistance*4;
%L(115).Resistance = L(115).Resistance*0.3;
%L(116).Resistance = L(116).Resistance*0.3;
%L(117).Resistance = L(117).Resistance*0.3;
%L(118).Resistance = L(118).Resistance*4;
%L(119).Resistance = L(119).Resistance*4;

%L([118 119]) = [];

%S(83).Resistance = 0.9;
%S(84).Resistance = 0.01;


%L([118 119]) = [];
% isub = find(strcmp({S.Name},'ELLERSLIE 89S'));
% indline1 = find(strcmp({L(:).fromSub},S(isub).Name));
% indline2 = find(strcmp({L(:).toSub},S(isub).Name));
% 
% indline = unique([indline1 indline2]);
% 
% for i = 1:length(indline)
%     if L(indline(i)).Voltage == 240
%         L(indline(i)).Resistance = 0.1*L(indline(i)).Resistance;
%     end
% end

%L(end+1) = L(117);
%L(end).Name = '909L #2';


%Remove HVDC Lines
L((find(contains({L(:).Name},'13L50')))) = [];
%L((find(contains({L(:).Name},'1325L')))) = [];
%L((find(contains({L(:).Name},'12L44')))) = [];
L((find(contains({L(:).Name},'12L41')))) = [];
%L((find(contains({L(:).Name},'1209')))) = [];
%L((find(contains({L(:).Name},'1202')))) = [];
%L((find(contains({L(:).Name},'1238')))) = [];
%L((find(contains({L(:).Name},'1239')))) = [];
%L((find(contains({L(:).Name},'12L70')))) = [];
%L((find(contains({L(:).Name},'12L85')))) = [];
%L((find(contains({L(:).Name},'1206')))) = [];

%L(234).Resistance = 0.5;
%S(153).Resistance = S(82).Resistance;
%T(152).W1 = 0.04;
%[L,S,T] = modify_network(L,S,T);

T(end+1).W1 = NaN; T(end).W2 = NaN; 
T(end).HV_Type = 'switch'; T(end).LV_Type = 'switch';
T(end).Sub = S(contains({S(:).Name},'LIVOCK 939S 500 kV')).Name;
T(end).W1Bus = 309;
T(end).W2Bus = 309;

%Remove May 11, 2024 Lines from Sherry Gao

% L((find(contains({L(:).Name},'928')))) = [];
% L((find(contains({L(:).Name},'938')))) = [];
% L((find(contains({L(:).Name},'9L07')))) = [];

% S((find(contains({S(:).Name},'ROSSDALE')))).Resistance = 10^9;
% S((find(contains({S(:).Name},'PETROLIA')))).Resistance = 10^9;
% S((find(contains({S(:).Name},'PAINTEARTH')))).Resistance = 10^9;
% S((find(contains({S(:).Name},'HORIZON MINING')))).Resistance = 10^9;



[L,S,T] = clean_network(L,S,T);

% toBus = L(contains({L(:).Name},'1202L')).toBus;
% fromBus = L(contains({L(:).Name},'1202L')).fromBus;
% 
% toSub = L(contains({L(:).Name},'1202L')).toSub;
% fromSub = L(contains({L(:).Name},'1202L')).fromSub;
% 
% L(contains({L(:).Name},'1202L')).toBus = fromBus;
% L(contains({L(:).Name},'1202L')).fromBus = toBus;
% 
% L(contains({L(:).Name},'1202L')).toSub = fromSub;
% L(contains({L(:).Name},'1202L')).fromSub = toSub;

%L(219).toBus = 309;
%L(220).fromBus = 309;

% for i = 1:length(S)
%     S(i).Resistance = 0.1;
% end

%T(154) = [];

% for i = 1:length(L)
%     L(i).Resistance = L(i).ResKm*L(i).Length/1000;
% end

%Perfect earthing case
% for i = 1:length(S)
%     if ~isnan(S(i).Resistance)
%         S(i).Resistance=10^-9; 
%     end
% end

% for i = 1:length(L)
%     L(i).Resistance = 0.1;
% end

%Turn off 240 kV lines
% for i = 1:nLines
%     if L(i).Voltage < 500
%         L(i).Resistance = 10^10;
%     end
% end

%Ground the Delta at the edge of the networks so that current can flow to
%ground at 240 kV buses. Does not change anything!
% for i = 1:nTrans
%     if isnan(T(i).W2)
%         T(i).W2 = 0.1;
%     end
% 
% end

%uniform = 1; latq= 0; lonq = 0; ex = 0; ey = 0; tind = 1;
[S,L,T,GIC_Subs,GIC_Trans,GIC_Lines,subLoc,nLines,nSubs,nTrans] = calc_gic_main(S,L,T,ex,ey,latq,lonq,tind,uniform);


%[Snew,Lnew,Tnew,g,t,l,subLoc2,nLines2,nSubs2,nTrans2] = calc_gic_main(Snew,Lnew,Tnew,ex,ey,latq,lonq,tind,uniform);
%%

disp(['***********COMPLETE: TOTAL TIME = ',num2str(toc(fullTic)),' seconds**************'])










