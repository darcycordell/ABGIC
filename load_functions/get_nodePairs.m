function [nodePairs,nodeRes,nodeInd,edges,indices,neutralNodes,autoind] = get_nodePairs(L,T,S)
%Get number of buses in the network by finding max bus number from list of
%buses
nBus = max([[T.W1Bus] [T.W2Bus] [L.fromBus] [L.toBus]]);

nNeutral = length(S); %default is 1 grounded neutral per substation
nNode = nBus + nNeutral; %total nodes in the network
nLines = length(L);
nTrans = length(T);

%Get substation indices for transformers
[~,tLoc] = ismember({T(:).Sub},{S(:).Name});

%Neutral node indices
neutralNodes = tLoc'+nBus;

%Pair up autotransformer nodes so that HV bus is connected
%directly to LV bus along series winding
autoind = find(strcmp({T(:).HV_Type},'auto'));
autoPairs = [[T(autoind).W1Bus]; [T(autoind).W2Bus]]';


%The autotransformer has confused me at various times. My understanding is
%that the transformer is wired together as this:
%
%
%       Node #1       tap       Series (W1)       Node #2
%  LV bus o------------.--------{}{}{}{}{}--------o HV bus
%                      |
%                      {}
%                      {}
%                      {} Common (W2)
%                      {}
%                      |
%                      |               ___
%            Node #3   o--------------- _   Ground point
%                                       .
%
% In this case:
%   The path from N1 to N3 only passes through W2 common resistance
%   The path from N2 to N1 only passes through W1 series resistance
%   The path from N2 to N3 passes through W1 and W2 resistances in series
%
%
tResHVAuto = [T(:).W1]; %resistance of HV windings
%autotransformers have path from HV bus to ground with series (W1) and
%common (W2) in series so resistances should add for auto transformers:
tResHVAuto(autoind) = [T(autoind).W1] + [T(autoind).W2]; 
%
% However, this method does not match Horton! In order for my results to
% match Horton, I have to make the path from N2 to N3 ignore the W1
% resistance and only take the W2 resistance:
tResHVAuto(autoind) = NaN;

% I don't know why this is the case since it doesn't make sense from any
% circuit diagram of an autotransformer I've ever seen, but I'm just
% leaving it as is for now since it matches Horton.

%Now I create a list of all nodes pair connections between buses
% First column is "from bus" and second column is "to bus"
% The first block is transmission line bus-to-bus pairs
% The second block is HV (series) bus node to neutral node
% The third block is LV (common) bus node to neutral node
% The fourth block is LV bus node to HV bus node for autotransformers
nodePairs = [[L(:).fromBus]' [L(:).toBus]'; ...
            [T(:).W1Bus]' neutralNodes; ...
            [T(:).W2Bus]' neutralNodes; ...
            autoPairs(:,1) autoPairs(:,2)];

%Corresponding resistances for each node pair:
nodeRes = [[L(:).Resistance] ... %Line resistances
            tResHVAuto ... % HV W1 resistances for wye, and NaN for W1 series in auto
            [T(:).W2] ... % LV W2 resistances for both wye and auto transformers
            [T(autoind).W1]]';  % W1 resistance for bus-to-bus auto transformers

%Corresponding indices for each node pair
nodeInd = [1:nLines 1:nTrans 1:nTrans autoind]';

%There is a way to make this faster by removing null rows prior to entering
%the time for loop, but it would take a fair amount of re-indexing all the
%matrices.
% indnull = isnan(nodeRes);
% nodePairs(indnull,:) = [];
% nodeRes(indnull) = [];
% nodeInd(indnull) = [];


%Get unique node pairs
[edges, ~, ia] = unique(nodePairs, 'rows','stable');

%Get indices of where each unique node pair appears in the original nodePairs variable
indices     = accumarray(ia, (1:numel(ia)).', [], @(r){sort(r)});