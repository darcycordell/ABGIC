function L = calc_line_resistance(L)

nLines = length(L);
for i = 1:nLines

    if isnan(L(i).ResKm)

    %If there is no data present for line resistance then:
    %Option #1:
        %I am using 0.0455 Ω/mile for 240 kV lines (same as Horton et al.
        % uses for 345 kV lines). Convert to 0.02827 Ω/km
        %I am using 0.0227 Ω/mile for 500 kV lines (same as Horton et al.
        % uses for 500 kV lines). Convert to 0.01411 Ω/km

%         if L(i).Voltage == min([L(:).Voltage])
%             L(i).Resistance = 0.02827*L(i).Length/1000;
%         elseif L(i).Voltage == max([L(:).Voltage])
%             L(i).Resistance = 0.01411*L(i).Length/1000;
%         else
%             error('Your power network has more than two voltage levels!')
%         end

    %Option #2: Use average value of lines of same voltage
        if L(i).Voltage == min([L(:).Voltage])
            ind = find([L(:).Voltage] == min([L(:).Voltage]));
            val = mean([L(ind).ResKm],'omitnan');

            L(i).ResKm = val;

            L(i).Resistance = val*L(i).Length/1000;

        elseif L(i).Voltage == max([L(:).Voltage])
            ind = find([L(:).Voltage] == max([L(:).Voltage]));
            val = mean([L(ind).ResKm],'omitnan');

            L(i).ResKm = val;

            L(i).Resistance = val*L(i).Length/1000;

        else
            error('Your power network has more than two voltage levels!')
        end

    else

        L(i).Resistance = L(i).ResKm*L(i).Length/1000;

    end

end