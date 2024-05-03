function L = get_downsampled_line(L,segLength)
%%
for i = 1:length(L)
    
    %Get individual transmission line path
    path = L(i).Loc;

    %Get line segments of original transmission line
    step = [0; distance(path(1:end-1,1),path(1:end-1,2),path(2:end,1),path(2:end,2),referenceEllipsoid('earth'))];
    
    %Re-sample each segment into multiple segLength segments.
    lineseglength = cumsum(step);
    
    %If you only do linspace with the start and end points, then you end
    %up with an empty vector if lineseglength < segLength
    Vq = unique([0 linspace(0,lineseglength(end), floor(lineseglength(end)/segLength)) lineseglength(end)]);
    
    %Re-sample each segment into 100 equally spaced segments.
    %Vq = linspace(0,lineseglength(end),100);
    
    
    interpline = interp1(lineseglength, path, Vq);
    
    %New line with vertices spaced every segLength meters.
    L(i).Loc = interpline;


end