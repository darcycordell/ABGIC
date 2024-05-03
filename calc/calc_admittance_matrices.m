function [Yn,Ye] = calc_admittance_matrices(edges,indices,nodeRes,neutralNodes,S)


Yn = zeros(max(edges(:))); %admittance matrix

for k = 1:size(edges,1)

    %All transformers connected to a bus are wired in parallel, so the
    %total resistance is summation 1/R_t = 1/R1 + 1/R2 + ...
    %and total admittance is Y_t = 1/R_t

    %Same is also true of double circuit lines in parallel

    nodalY = nansum(-1./nodeRes(indices{k})); % = 1/R_t
    
    Yn(edges(k,1),edges(k,2)) = nodalY;
    Yn(edges(k,2),edges(k,1)) = nodalY;

end

%Get diagonal elements as sum of all other elements in Yn
diagYn = -sum(Yn,2);

%Add the diagonal elements to the Yn matrix
Yn = Yn + diag(diagYn);

%Construct earthing admittance matrix
diagYe = zeros(length(diagYn),1);
diagYe(unique(neutralNodes))=(1./[S(~isnan([S(:).Resistance])).Resistance])/3; %divide *admittance* by 3 for 3 phase power
                                       %This is the same as multiplying the
                                       %resistance by 3 (pers. comm. David
                                       %Boteler, Jan 4, 2024)

Ye = diag(diagYe);

%Option to include effects of nearby substations to earthing matrix
%Does not make any noticeable difference
% subLoc = reshape([S(:).Loc],2,length(S))';
% Lij = [];
% for i = 1:length(S)
%     for j = 1:length(S)
%         Lij(i,j) = distance(subLoc(i,1),subLoc(i,2),subLoc(j,1),subLoc(j,2),referenceEllipsoid('WGS84'));
%     end
% end
% 
% Ls = min(Lij(Lij~=0));
% 
% nY = size(Ye,1);
% k = nY-length(S);
% for i = 1:length(S)
%     for j = 1:length(S)
% 
%         if i~=j
% 
%             Z = ((1/Ye(i+k,i+k)+1/Ye(j+k,j+k))/4)*(Ls/Lij(i,j));
%             Ye(i+k,j+k) = Z;
%         end
%     end
% end
% 
% 
