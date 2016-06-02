function dU = grad_potential(b,x0,xT,obs_noise,f,noise_level,df,dt)
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
   x(i) = x(i-1)+dt*f(x(i-1)) + b(i-1); 
end

%the gradient of the path variables
dX = zeros(n-1,1);

dX(n-1) = 1;

for i=2:n-1
   dX(n-i) = dX(n-i+1) * (1+df(x(n-i+1))*dt);
end

dU = -(xT-x(n))/obs_noise*dX + b/(noise_level^2*dt);