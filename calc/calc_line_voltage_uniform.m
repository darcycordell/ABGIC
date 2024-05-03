function V = calc_line_voltage_uniform(ex,ey,L,S)
% Get Line Voltages-------------------------------------------------------


[~,sub1] = ismember({L(:).fromSub},{S(:).Name});
[~,sub2] = ismember({L(:).toSub},{S(:).Name});

subLoc = reshape([S(:).Loc],2,length(S))';

phi = (subLoc(sub1,1)+subLoc(sub2,1))/2; %average latitude


[~,az] = distance(subLoc(sub1,1),subLoc(sub1,2),subLoc(sub2,1),subLoc(sub2,2),referenceEllipsoid('WGS84'));
[Dn] = distance(subLoc(sub1,1),subLoc(sub1,2),subLoc(sub2,1),subLoc(sub1,2),referenceEllipsoid('WGS84'));
[De] = distance(phi,subLoc(sub1,2),phi,subLoc(sub2,2),referenceEllipsoid('WGS84'));

%azind1 = all([az>=0 az<90],2);
azind2 = all([az>=90 az<180],2); 
azind3 = all([az>=180 az<270],2);
azind4 = all([az>=270 az<360],2);

Dn(azind2) = -Dn(azind2);
Dn(azind3) = -Dn(azind3);
De(azind3) = -De(azind3);
De(azind4) = -De(azind4);

%Line voltages
V = zeros(length(ex),length(Dn));
if length(ex) >= length(Dn)

    for i = 1:length(Dn)
        V(:,i) = (ex*Dn(i)+ey*De(i))/1000;
    end

else

    for i = 1:length(ex)
        V(i,:) = (ex(i)*Dn+ey(i)*De)/1000;
    end


end