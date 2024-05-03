function [GIC_Subs,GIC_Lines,GIC_Trans] = calc_gic(L,T,V,Yn,Ye,nodePairs,nodeRes,autoind,indices,edges,indnull,indnotnull,nBus)

nLines = length(L);
nTrans = length(T);

%Current sources in network edges are zero everywhere
J = zeros(size(nodeRes));
%except in transmission lines:
J(1:nLines) = V./[L(:).Resistance];

%Driving emf in network edges is zero everywhere
emf = zeros(size(nodeRes));
%except in transmission lines:
emf(1:nLines) = V;

nNode = size(Yn,1);
Je = zeros(nNode); %current matrix
Vm = zeros(nNode); %voltage matrix
for k = 1:size(edges,1)

    %Get total current flowing in/out of nodes
    sumcurrent = sum(-J(indices{k}),'omitnan');
    
    %Fill in current matrix
    Je(edges(k,1),edges(k,2)) = sumcurrent;
    Je(edges(k,2),edges(k,1)) = -sumcurrent;

    %And voltage matrix
    Vm(edges(k,1),edges(k,2)) = mean(emf(indices{k}),'omitnan');
    Vm(edges(k,2),edges(k,1)) = -mean(emf(indices{k}),'omitnan');


end

%Remove the nodes that do not have GIC flow (e.g. termination points, delta
%windings, switching stations, ungrounded nodes)
Yn(indnull,:) = [];
Yn(:,indnull) = [];

Ye(indnull,:) = [];
Ye(:,indnull) = [];

Je(indnull,:) = [];
Je(:,indnull) = [];

Vm(indnull,:) = [];
Vm(:,indnull) = [];

Js = sum(Je,2); %Nodal current sources

%Calculate nodal voltages
Vn = (Yn+Ye)\Js; %Equation 22 of Pirjola et al. 2021.

%Calculate nodal voltages plus driving EMF
Vmat = Vm + (Vn-Vn');


Ie = Ye*Vn; %grounding currents (per phase), Equation 18 of Pirjola et al. 2021
            %Need to multiply by three to get the true flow to ground

groundRes = zeros(nNode,1);
groundRes(indnotnull) = Ie;


%Ikn = Yn.*Vmat; %compute current flowing in/out of nodes (Equation 2 of Pirjola et al, 2021)

[Vmat] = insert_arbitrary_RowsColumns(Vmat,indnull);

%Now its just a matter of finding the Ikn that match with transformers and
%transmission lines
Ival = diag(Vmat(nodePairs(:,1),nodePairs(:,2))./nodeRes);


GIC_Subs = groundRes(nBus+1:end);
GIC_Lines = Ival(1:nLines);
GIC_Trans_W1 = Ival(nLines+1:nLines+nTrans);
GIC_Trans_W2 = Ival(nLines+nTrans+1:nLines+2*nTrans);
GIC_Trans_W1(autoind) = Ival(nLines+2*nTrans+1:end);

GIC_Trans = [GIC_Trans_W1 GIC_Trans_W2];


