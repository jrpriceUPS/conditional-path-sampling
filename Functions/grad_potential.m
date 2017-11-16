function dU = grad_potential(path,start_dist,end_dist,f,df,noise_level,dt)
%Takes the path and computes the potential dU
%
%INPUTS:
%
%path         =  the path
%start_dist   =  the negative logarithm of the starting distribution
%end_dist     =  the negative logarithm of the ending distribution
%f            =  the drift function
%f            =  the derivative of the drift function
%noise_level  =  the level of Brownian noise (sigma)
%dt           =  the time increment

%compute the length of the path
n = length(path);

%initialize the gradient vector (see May note)
dU = zeros(n,1);

%compute the gradient of the first point (see May note)
dU(1) = -(path(2)-(path(1)+f(path(1))*dt))/(noise_level^2*dt).*(1+df(path(1))*dt) + start_dist(path(1));

%compute the intermediate gradients (see May note)
dU(2:n-1) = (path(2:n-1)-(path(1:n-2)+f(path(1:n-2))*dt))/(noise_level^2*dt)...
            - (path(3:n)-(path(2:n-1)+f(path(2:n-1))*dt))/(noise_level^2*dt).*(1+df(path(2:n-1))*dt);
      
%compute the ending point gradient (see May note)
dU(n) = (path(n) - (path(n-1)+f(path(n-1))*dt))/(noise_level^2*dt) + end_dist(path(end));