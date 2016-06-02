function U = potential(path,start_dist,end_dist,f,noise_level,dt)
%Takes the path and computes the potential U under the explicit
%Euler-Maruyama scheme
%
%INPUTS:
%
%path         =  the path
%start_dist   =  the negative logarithm of the starting distribution
%end_dist     =  the negative logarithm of the ending distribution
%f            =  the drift function
%noise_level  =  the level of Brownian noise (sigma)
%dt           =  the time increment

%compute the potential (see May note)
U = sum((path(2:end)-(path(1:end-1)+f(path(1:end-1))*dt)).^2/(2*noise_level^2*dt))...
    + start_dist(path(1))...
    + end_dist(path(end));