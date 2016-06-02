function U = periodic_potential_v2(path,x0,xT,obs_noise,f,noise_level,dt,L,R)
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

%compute distance between end point and conditional ending mean
distance = min([abs(xT-path(end)),abs(xT+(R-L)-path(end)),abs(xT-(R-L)-path(end))]);

%compute the potential
U = sum((path(2:end)-(path(1:end-1)+f(path(1:end-1))*dt)).^2/(2*noise_level^2*dt))...
    +distance^2/(2*obs_noise^2);