function V = calc_line_voltage(L,latq,lonq,ex,ey,method)
%Calculate line integral along transmission line paths given vector data
%
% Inputs: L is the line network structure. Each entry of L includes the
%           name, buses, transmission voltage, length, line resistance, and
%           line coordinates. Line coordinates have already been sorted
%           from the fromBus to the toBus
%         latq, lonq: the coordinates where the vector data are defined
%         ex, ey: Vector data in x (north) and y (east). Electric fields
%               are in V/m
%         method is the interpolation method (e.g. cubic, nearest neighbor)
%
%
% Outputs:
%       V is the (ntimes x nlines) matrix of line voltages in Volts
%
%There are two different methods to do the interpolation of the E-fields
%onto the transmission lines. The first uses a FEX script called func2mat
%which does the interpolation of *all* transmission line points and all
%time steps in a single step. This is *incredibly* fast for large problems
%and allows for the calculation of the line voltages for a 10000s
%of points in tens of seconds.
%
%The second option does the interpolation in nested for loops, looping over
%each transmission line and each time interval. This is incredibly slow
%when considering many time steps. For example, 1000s of time steps would
%take 10s of *minutes*. However, there is a tradeoff somewhere, because
%this method is much faster if you are only doing a few timesteps. For
%example, doing only 1 timestep takes only 2 seconds, but still takes about
%35 seconds using the other method. 
%
%I tested this a bit and there doesn't seem to be any hard and fast rule
%for which is faster. But, roughly it seems that if you are doing >500 time
%steps, then the first option is faster, so I've put an if statement to
%check this.
%
% Note that the two methods give identical results within numerical error.

nTimes = size(ex,1);
nLines = length(L);

if nTimes>=500
    speedup_option = 1;
else
    speedup_option = 0;
end

speedup_option = 1;

%Initialize large matrices
V = zeros(nTimes,nLines);


%The origin of the survey used for lat-long to UTM conversions
% origlon =(max(lonq)-min(lonq))/2+min(lonq);
% origlat = (max(latq)-min(latq))/2+min(latq);

rx = []; ry = []; indLines = {}; midlat = []; midlon = [];
for sidx = 1:nLines
    
    indLines{sidx} = 1+length(rx):length(rx)+length(L(sidx).Loc(:,1));
    
    rx = [rx L(sidx).Loc(:,1)'];
    ry = [ry L(sidx).Loc(:,2)'];

    np = length(L(sidx).Loc(:,1));

    %Find midpoint of each line for later UTM conversion
    midlat = [midlat L(sidx).Loc(round(np/2),1)*ones(1,np)];
    midlon = [midlon L(sidx).Loc(round(np/2),2)*ones(1,np)];
   
end

disp('Finished collecting all transmission line coordinates')

 %Two ways to do the conversion to UTM. First way uses the center point of
 %the entire survey area, whereas the second way uses the midpoint of each
 %line segment. Using the midpoint leads to less error in the conversion
 %due to projection issues. The difference between the two can actually be
 %quite significant (e.g. up to 4 V)
%[y,x] = geo2utm(ry,rx,origlon,origlat); %convert using the full dataset midpoint as an origin (x = northing, y = easting)
[y,x] = geo2utm(ry,rx,midlon,midlat); %convert using the full dataset midpoint as an origin (x = northing, y = easting)
  
Fex = scatteredInterpolant(latq,lonq,ex(1,:).',method);
Fey = scatteredInterpolant(latq,lonq,ey(1,:).',method);  

if speedup_option == 1 %"faster" option

    disp('nTimes>500. Using func2mat method Y(:)=A*X(:)----')
    
    disp('...interpolating E fields onto transmission line coordinates........')
    disp('Working on Ex........')
    Ax=func2mat(@(D) func_func2mat(D,Fex,rx,ry), ex(1,:).');
    disp('Finished Ex........')

    disp('Working on Ey........')
    Ay=func2mat(@(D) func_func2mat(D,Fey,rx,ry), ey(1,:).');
    disp('Finished Ey........')
    
    
    su_window = 3600; %to save memory, I do a loop through a set of time-steps (e.g. 1 hour)
            %I can do the matrix multiplication in one step, but this
            %results in exnn3d, eynn3d, exnn1d, and eynn1d matrices that
            %are huge and exceed 16 GB of RAM. I can save memory by doing
            %each hour at a time
    su_ind = unique([0:su_window:nTimes nTimes]);

    
    for su = 1:length(su_ind)-1
        
        ex_int=Ax*ex(su_ind(su)+1:su_ind(su+1),:).';
        ey_int=Ay*ey(su_ind(su)+1:su_ind(su+1),:).';
        
        for sidx = 1:nLines

            s = indLines{sidx};

            indslice = su_ind(su)+1:su_ind(su+1);

            dx = diff(x(s));
            dy = diff(y(s));

            d = distance(rx(s(1)),ry(s(1)),rx(s(end)),ry(s(end)),referenceEllipsoid('WGS84'));

            dex = (ex_int(s(1:end-1),:)+ex_int(s(2:end),:))/2;
            dey = (ey_int(s(1:end-1),:)+ey_int(s(2:end),:))/2;

            if length(s)>1
                V(indslice,sidx) = nansum(repmat(dx,size(dex,2),1).*dex.'+repmat(dy,size(dey,2),1).*dey.',2);
            end

        end
        disp([num2str(100*su_ind(su+1)/nTimes,'%.1f'),'% of Time Steps Completed.............................'])

    end


else %Original option using nested for loops

    disp('nTimes<500. Looping over each time step-----')

    V = zeros(nTimes,nLines);

    for sidx = 1:nLines %Loop through lines

        

        %Line vertices in latitude (rx) and longitude (ry)
        rx = L(sidx).Loc(:,1)';
        ry = L(sidx).Loc(:,2)';

        %Convert lat long to metres (UTM). Note: y = easting; x = northing 
        % Similar to above, there's two ways to do it, either using the
        % midpoint of the line, or using the center point of the survey. I
        % prefer to use the midpoint since it leads to less projection
        % error.
        [y,x] = geo2utm(ry,rx,ry(round(length(ry)/2)),rx(round(length(rx)/2))); %convert using the transmission line segment midpoint as a reference
        %[y,x] = geo2utm(ry,rx,origlon,origlat); %convert using the full dataset midpoint as an origin (x = northing, y = easting)
        tic

        %Interpolate transmission line values and convert to V/m
        % Here I use scatteredInterpolant to generate 4 functions for ex1d,
        % ey1d, ex3d, and ey3d. Once I have the interpolant function, then I
        % enter the loop and just replace the function values since d.loc never
        % changes. This is faster than e.g. gridding on every time step.

        Fex = scatteredInterpolant(latq,lonq,ex(1,:).',method);
        Fey = scatteredInterpolant(latq,lonq,ey(1,:).',method);    


        for tidx = 1:nTimes %Loop through each time step

            Fex.Values = ex(tidx,:).';
            ex_int = Fex(rx,ry);

            Fey.Values = ey(tidx,:).';
            ey_int = Fey(rx,ry);


    %         Perform line integral (x and y are in meters, e is in V/m)
    %         Resources:
    %               https://www.mathworks.com/matlabcentral/answers/441416-numerical-calculation-of-line-integral-over-a-vector-field
    %               https://ocw.mit.edu/ans7870/18/18.013a/textbook/HTML/chapter25/section04.html
            V(tidx,sidx) = (nansum((diff(x)).*(ex_int(1:end-1)+ex_int(2:end))/2+(diff(y)).*(ey_int(1:end-1)+ey_int(2:end))/2));


        end

        toc
        disp(['Transmission Line #',num2str(sidx),' of ',num2str(nLines),' Completed...........................'])

    end

end
