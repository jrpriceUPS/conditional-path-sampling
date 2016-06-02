function dU = periodic_grad_potential_v2(path,x0,xT,obs_noise,f,noise_level,df,dt,L,R)
%Takes the sequence of increments and computes the potential U and the 
%gradient of the potential dU for a one dimensional problem
%
%INPUTS:
%
%b      =  the sequence of increments
%
%x0     =  the initial condition
%
%xT     =  the mean of the end gaussian
%
%sigma  =  the variance of the end gaussian
%
%f      =  the drift function
%
%df     =  the derivative of the drift functoin
%
%dt     =  the time increment

%append the initial condition
path = [x0;path];

n = length(path);

dU = zeros(n-1,1);

distance_1 = zeros(n-2,3);

distance_1(:,1)   =  path(2:n-1)-(path(1:n-2)+f(path(1:n-2))*dt);
distance_1(:,2)   =  -distance_1(:,1) + (R-L);
distance_1(:,3)  =  distance_1(:,1) + (R-L);

[~,idx1]=min(abs(distance_1),[],2);

mask1 = [idx1==1 idx1==2 idx1==3];
d1 = sum(distance_1.*mask1,2);


distance_2 = zeros(n-2,3);

distance_2(:,1)   =  path(3:n)-(path(2:n-1)+f(path(2:n-1))*dt);
distance_2(:,2)   =  distance_2(:,1) - (R-L);
distance_2(:,3)   =  distance_2(:,1) + (R-L);

[~,idx2]=min(abs(distance_2),[],2);

mask2 = [idx2==1 idx2==2 idx2==3];
d2 = sum(distance_2.*mask2,2);


distance_3 = zeros(1,3);

distance_3(1)   =  path(n) - (path(n-1)+f(path(n-1))*dt);
distance_3(2)   =  distance_3(1) - (R-L);
distance_3(3)   =  distance_3(1) + (R-L);

[~,idx3]=min(abs(distance_3));

d3 = distance_3(idx3);


%compute distance between end point and conditional ending mean
end_distance = zeros(1,3);

end_distance(1) = xT-path(n);
end_distance(2) = xT-path(n) - (R-L);
end_distance(3) = xT-path(n) + (R-L);

[~,idx_end]=min(end_distance);

final = end_distance(idx_end);


dU(1:n-2) = d1/(noise_level^2*dt)...
            - d2/(noise_level^2*dt).*(1+df(path(2:n-1))*dt);
        
dU(n-1) = d3/(noise_level^2*dt) - final/(obs_noise^2);