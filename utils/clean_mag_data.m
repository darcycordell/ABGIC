function [xclean,yclean,count] = clean_mag_data(x,y,fs)

% %Remove "bad" spikes where the derivative is >0.05
% indx = find((Bx_all(1:end-1)-Bx_all(2:end))./Bx_all(1:end-1)>0.05);
% indy = find((By_all(1:end-1)-By_all(2:end))./By_all(1:end-1)>0.05);
% ind = unique([indx; indy]);
% 
% %Meanook site has a large chunk of "bad" data from 15:09 to 15:14 UT which
% %is not a "spike" so is missed by the above derivative method.
% if strcmp(site,'MEA') 
%     ind = [ind; [54580:54860]'-tstart];
% end



count = 0;
while 1
%Remove "bad" spikes where the derivative is >0.05
dbxdt = abs((x(1:end-1)-x(2:end))./(1/fs));
dbydt = abs((y(1:end-1)-y(2:end))./(1/fs));
indx = find(dbxdt>500);
indy = find(dbydt>500);
ind = unique([indx; indy]);

count = count + length(ind);

x(ind) = NaN; 
y(ind) = NaN; 

if isempty(ind)
    break
end
% 
% %Meanook site has a large chunk of "bad" data from 15:09 to 15:14 UT which
% %is not a "spike" so is missed by the above derivative method.
% if strcmp(site,'MEA')
%     tind = find(times_all==datetime('2017-09-08 15:08:00')):find(times_all==datetime('2017-09-08 15:15:00')); 
%     ind = [ind; tind'];
% end

end

%Interpolate NaNs
xclean = inpaint_nans(x);
yclean = inpaint_nans(y);



