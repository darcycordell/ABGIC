function [L,T] = get_buses(L,S)

nSubs = length(S);
nBus = nSubs*2; %assume two buses per substation

indnan = find(isnan([L.Voltage]));
L(indnan) = [];
num_voltage_levels = length(unique([L.Voltage]));

if num_voltage_levels ~= 2
    error('Your network does not have 2 voltage levels...')
end

LV = min(unique([L.Voltage]));
HV = max(unique([L.Voltage]));

L(1).('fromBus') = [];
L(1).('toBus') = [];

count = 1;
for i = 1:nSubs

    indLinesFrom = find(strcmp(S(i).Name,{L(:).fromSub}));
    indLinesTo = find(strcmp(S(i).Name,{L(:).toSub}));
    voltage_levels = [L(indLinesFrom).Voltage L(indLinesTo).Voltage];

    if length(unique(voltage_levels)) == 1 %HV has lines, LV is termination

        T(i,1).W1 = 0.1;
        T(i,1).W2 = NaN;
        T(i,1).HV_Type = 'wye';
        T(i,1).LV_Type = 'delta';
        T(i,1).Sub = S(i).Name;
        
        T(i,1).W1Bus = count;
        T(i,1).W2Bus = count+1;
    

        for j = 1:length(indLinesFrom)
            L(indLinesFrom(j)).fromBus = T(i).W1Bus;
        end

        for j = 1:length(indLinesTo)
            L(indLinesTo(j)).toBus = T(i).W1Bus;
        end

        

    else %both HV and LV buses have lines, step-down transformer

        T(i,1).W1 = 0.04; %series
        T(i,1).W2 = 0.06; %common
        T(i,1).HV_Type = 'auto';
        T(i,1).LV_Type = 'auto';
        T(i,1).Sub = S(i).Name;
        
        T(i,1).W1Bus = count;
        T(i,1).W2Bus = count+1;

        for j = 1:length(indLinesFrom)
            if L(indLinesFrom(j)).Voltage == HV
                L(indLinesFrom(j)).fromBus = T(i).W1Bus;
            else
                L(indLinesFrom(j)).fromBus = T(i).W2Bus;
            end
        end

        for j = 1:length(indLinesTo)
            if L(indLinesTo(j)).Voltage == HV
                L(indLinesTo(j)).toBus = T(i).W1Bus;
            else
                L(indLinesTo(j)).toBus = T(i).W2Bus;
            end
        end

    end

    
    count = count+2;


end




   

