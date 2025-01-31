function [S,L,T,GIC_Subs,GIC_Trans,GIC_Lines,subLoc,nLines,nSubs,nTrans,V] = calc_gic_main(S,L,T,ex,ey,latq,lonq,tind,uniform)
subLoc = reshape([S(:).Loc],2,length(S))';
nLines = length(L);
nSubs = length(S);
nTrans = length(T);


%Compute E field at each substation
% 
% [bxSub,bySub] = interpolate_b(b,subLoc(:,1),subLoc(:,2));
% BxSub = zeros(length(latq),nf); BxSub = BxSub+1i*BxSub;
% BySub = BxSub;
% for i = 1:size(bxSub,1)
%     [BxSub(i,:),BySub(i,:)] = calc_fft(bxSub(i,:)',bySub(i,:)',b(1).fs,b(1).pad); 
% end
% 
% latQuery = subLoc(:,1); lonQuery = subLoc(:,2);
% [~,indMT] = min(distance(latQuery-zeros(1,d.ns),lonQuery-zeros(1,d.ns),d.loc(:,1)'-zeros(nSubs,1),d.loc(:,2)'-zeros(nSubs,1),referenceEllipsoid('wgs84')),[],2); 
% [exSub,eySub] = get_e(d.Z(:,:,indMT),d.f,BxSub,BySub,b);
% 


%thetaVec = [0:5:360];
thetaVec = 0;
%thetaVec = 117.3227;



disp('************************LINE VOLTAGES****************************')

if ~uniform
    V = calc_line_voltage(L,latq,lonq,ex(tind,:),ey(tind,:),'natural');
    %V1 = calc_line_voltage(L,latq,lonq,ex1(tind,:),ey1(tind,:),'natural');
    %V(234) = -V(234);
    nTimes = size(V,1);

else
    nTimes = 1;

end


disp('********************LINE VOLTAGES:COMPLETED************************')



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
%s1d = GIC_Subs; l1d = GIC_Lines; t1d = GIC_Trans;
%su = GIC_Subs; lu = GIC_Lines; tu = GIC_Trans;

for thetaidx = 1:length(thetaVec)
    
if uniform
    
    %[su,lu,tu] = calc_gic(L,T,Vu,Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);
    x = cosd(thetaVec(thetaidx));
    y = sind(thetaVec(thetaidx));
    x = 0;
    y = 1;
    Vu = calc_line_voltage_uniform(x,y,L,S);
    [GIC_Subs(:,thetaidx),GIC_Lines(:,thetaidx),GIC_Trans(:,:,thetaidx)] = calc_gic(L,T,Vu,Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);
        

else
ticHour = tic; ticMin = tic;
for i = 1:nTimes

    %mageuni = 0.5118;
    %mageuni = 1;
    %x = mageuni*cosd(thetaVec(thetaidx));
    %y = mageuni*sind(thetaVec(thetaidx));
    %Vu = calc_line_voltage_uniform(x,y,L,S);
    %%
    %V = [222.93 0 -222.93];
    [GIC_Subs(:,i),GIC_Lines(:,i),GIC_Trans(:,:,i)] = calc_gic(L,T,V(i,:),Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);
    
    %[s1d(:,i),l1d(:,i),t1d(:,:,i)] = calc_gic(L,T,V1(i,:),Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);
    
    %[su(:,i),lu(:,i),tu(:,:,i)] = calc_gic(L,T,Vu(i,:),Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus);
    

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

% g(:,thetaidx) = su;
% l(:,thetaidx) = lu;
% t(:,:,thetaidx) = tu;
% v(:,thetaidx) = Vu./[L(:).Resistance];


end