function U = periodic_potential(b,x0,xT,obs_noise,f,sigma,dt,L,R)
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

%length of the path
n = length(b)+1;

%the path itself
x = zeros(n,1);
x(1) = x0;

%compute the path
for i=2:n
   x(i) = x(i-1)+dt*f(x(i-1))+b(i-1); 
   
   if x(i) < L
       x(i) = R-(L-x(i));
       
   elseif x(i) > R
       
       x(i) = x(i) - R + L;
   end
end

%compute the potential
U = (xT - x(n))^2/(2*obs_noise) + sum(b.^2/(2*sigma^2*dt));