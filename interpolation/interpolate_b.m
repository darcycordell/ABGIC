function [xint,yint] = interpolate_b(b,latq,lonq)


s_lat = [b(:).lat]; %latitude points of sites
s_lon = [b(:).lon]; %longitude points of sites
x = [b(:).x]; %bx data
y = [b(:).y]; %by data


% Interpolate onto query points

t_ind = 1;

%Create scattered interpolant using single time step
Fx_int = scatteredInterpolant(s_lon',s_lat',x(t_ind,:)','natural');
Fy_int = scatteredInterpolant(s_lon',s_lat',y(t_ind,:)','natural');

%tic
%Super fast way to compute interpolated field at all query points and all
%time steps in a single step without for loops:
Ax=func2mat(@(D) func_func2mat(D,Fx_int,lonq,latq), x(t_ind,:)' );
xint=Ax*x';

Ay=func2mat(@(D) func_func2mat(D,Fy_int,lonq,latq), y(t_ind,:)' );
yint=Ay*y';
%toc
