clearvars

get_files;

cd(mainPath);

fullTic = tic;



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
subLoc = reshape([S(:).Loc],2,length(S))';
nLines = length(L);
nSubs = length(S);
nTrans = length(T);

%Perfect earthing case
%for i = 1:length(S); S(i).Resistance=10^-9; end

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
%%


[d] = get_Z(1); %Option flag 1 = MT impedance
                %Option flag 2 = 1000 Ωm halfspace


[GIC_Hall,tgic] = get_hall_data;
tgic_orig = tgic;


%Compute E field at each impedance query point
latq = d.loc(:,1);
lonq = d.loc(:,2);


[bx,by,Bx,By,b] = get_b(latq,lonq);
[ex,ey] = get_e(d.Z,d.f,Bx,By,b);


%Compute E field at each substation
[bxSub,bySub,BxSub,BySub] = get_b(subLoc(:,1),subLoc(:,2));

latQuery = subLoc(:,1); lonQuery = subLoc(:,2);
[~,indMT] = min(distance(latQuery-zeros(1,d.ns),lonQuery-zeros(1,d.ns),d.loc(:,1)'-zeros(nSubs,1),d.loc(:,2)'-zeros(nSubs,1),referenceEllipsoid('wgs84')),[],2); 
[exSub,eySub] = get_e(d.Z(:,:,indMT),d.f,BxSub,BySub,b);


%t1 = find(b(1).times==datetime('2021-10-12 10:54:00'));
[maxe, maxeIndices] = max(sqrt(ex.^2+ey.^2));
[maxeCounts, maxeGroups] = groupcounts(maxeIndices');
tind = maxeGroups(find(maxeCounts==max(maxeCounts)));

t1 = datetime('2023-04-23 23:00:00');
t2 = datetime('2023-04-24 10:00:00');

%t1 = datetime('2021-11-04 06:00:00');
%t2 = datetime('2021-11-04 13:00:00');

%t1 = datetime('2023-03-23 04:00:00');
%t2 = datetime('2023-03-23 18:00:00');

% t1 = datetime('2023-04-24 2:00:00');
% t2 = datetime('2023-04-24 3:00:00');


tind = isbetween(b(1).times,t1,t2);
%tind = find(b(1).times=='2023-04-24 07:26:46');
%tind = 1:length(b(1).times);




%%
disp('************************LINE VOLTAGES****************************')
V = calc_line_voltage(L,latq,lonq,ex(tind,:),ey(tind,:),'natural');
%V = calc_line_voltage_uniform(1,0,L,S);

disp('********************LINE VOLTAGES:COMPLETED************************')


nTimes = size(V,1);
GIC_Times = b(1).times(tind);
GIC_Times_orig = GIC_Times;

%%
nBus = max([[T.W1Bus] [T.W2Bus] [L.fromBus] [L.toBus]]);

%End Lines and Subs--------------------------------------------------------

disp('**************************NETWORK TOPOLOGY*************************')
%This function sets up the model topology
disp('...setting up nodes')
[nodePairs,nodeRes,nodeInd,edges,indices,neutralNodes,autoind] = get_nodePairs(L,T,S);

%This function gets admittance matrices
disp('...building admittance matrices')
[Yn,Ye] = calc_admittance_matrices(edges,indices,nodeRes,neutralNodes,S);

%Find indices of nodes that do not have GIC flow
indnull = find(diag(Yn)==0);
indnotnull = find(diag(Yn)~=0);


%------------------DONE MODEL SETUP--------------------------------------
% Everything above can be done without time loop


disp('***********************CALCULATE GIC*************************')
%Final step: Compute GIC in network:
tic
GIC_Subs = zeros(nSubs,nTimes); %Amps *per phase* flowing to ground
GIC_Lines = zeros(nLines,nTimes); %Amps *per phase* in lines
GIC_Trans = zeros(nTrans,2,nTimes); %Amps *per phase* in each winding
ticHour = tic; ticMin = tic;
parfor i = 1:size(V,1)
    [GIC_Subs(:,i),GIC_Lines(:,i),GIC_Trans(:,:,i)] = calc_gic(L,T,V(i,:),Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);

%     if mod(i,3600) == 0
%         disp([num2str(i/3600),' Hours Completed in ',num2str(toc(ticHour)),' seconds...............................'])
%         ticHour = tic;
%     end

%     if mod(i,60) == 0
%         disp([num2str(i/60),' Minutes Completed in ',num2str(toc(ticMin)),' seconds'])
%         ticMin = tic;
%         
%     end

i

end
toc

[maxSub,indmaxSub] = max(abs(GIC_Subs),[],2);
[~,indsortmaxSub] = sort(maxSub,'descend');
indsortmaxSub(isnan(maxSub(indsortmaxSub)))=[];


for i = 1:nLines
    L(i).GIC = GIC_Lines(i,:);
end

for i = 1:nSubs
    S(i).GIC = GIC_Subs(i,:);

    S(i).maxGIC = maxSub(i);
   
end

disp('*******************COMPLETED: GIC***************************')
disp('')
disp('')

fprintf(['Node #1       tap       Series (W1)       Node #2\n' ...
  'LV bus o------------.--------{}{}{}{}{}--------o HV bus\n' ...
  '                    |\n' ...
  '                    {}\n' ...
  '                    {}\n' ...
  '                    {}  Common (W2)\n' ...
  '                    {}\n' ...
  '                    |\n' ...
  '                    |               ___\n' ...           
  '          Node #3   o--------------- _   Ground point\n' ...
  '                                     . \n\n\n'])

disp(['***********COMPLETE: TOTAL TIME = ',num2str(toc(fullTic)),' seconds**************'])









