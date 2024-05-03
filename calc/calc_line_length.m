function L = calc_line_length(L,S)

nLines = length(L);

for i = 1:nLines

    if ~isempty(L(i).Loc)
        
        %Find indices of to and from substations
        indFrom = find(strcmp(L(i).fromSub,{S.Name}));
        indTo = find(strcmp(L(i).toSub,{S.Name}));
    
        %Get first line coordinate
        latFirst = L(i).Loc(1,1);
        lonFirst = L(i).Loc(1,2);
    
        %Compute the distance between the first line coordinate and the to/from
        %substations
        [d1] = distance(latFirst,lonFirst,S(indFrom).Loc(1),S(indFrom).Loc(2),referenceEllipsoid('WGS84'));
        [d2] = distance(latFirst,lonFirst,S(indTo).Loc(1),S(indTo).Loc(2),referenceEllipsoid('WGS84'));
    
        %If the line coordinates are going fromSub to toSub, then d1 < d2 and
        %we don't need to swtich the coordinates. However, if
        if d1 > d2
    
            %Then, this means the coordinates are going from the "toSub" index
            %to the "fromSub" index, so we need to flip the coordinates:
            L(i).Loc = flipud(L(i).Loc);
    
        end
    
        %Now, once the lines are oriented the right way, we can add the
        %substation locations as the first and last point of the line
        %coordinates:
    
        L(i).Loc = [S(indFrom).Loc; L(i).Loc; S(indTo).Loc];
    
    
        %Now once we have the full set of coordinates, we compute the distance
        %and resistance. Add 3% of the distance to account for sag.
    
        L(i).Length = 1.03*sum(distance(L(i).Loc(1:end-1,1),L(i).Loc(1:end-1,2), ...
            L(i).Loc(2:end,1),L(i).Loc(2:end,2),referenceEllipsoid('wgs84')));

    end
    
end