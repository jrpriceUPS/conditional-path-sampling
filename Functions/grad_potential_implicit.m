function dU = grad_potential_implicit(path,start_dist,end_dist,f,noise_level,df,df2,dt)
%Takes the path and computes the gradient of the potential U under the 
%implicit Euler-Maruyama scheme
%
%INPUTS:
%
%path         =  the path
%start_dist   =  the negative logarithm of the starting distribution
%end_dist     =  the negative logarithm of the ending distribution
%f            =  the drift function
%df            =  the derivative of the drift function
%df2          =  the second derivative of the drift function
%noise_level  =  the level of Brownian noise (sigma)
%dt           =  the time increment

%compute the length of the path
n = length(path);

%initialize the gradient vector
dU = zeros(n,1);

%note sliding windows of the path
x_0 = path(1:n-2);
x_1 = path(2:n-1);
x_2 = path(3:n);

%evaluate the drift and its derivatives for windows
f_0 = f(x_0);
f_1 = f(x_1);
df_0 = df(x_0);
df_1 = df(x_1);
df2_1 = df2(x_1);

%compute the gradient of the first point
dU(1)  =  -((path(2)-path(1))*(1-df_0(1)*dt)-f_0(1)*dt)/(noise_level^2*dt)*(1+(path(2)-path(1))*df2(path(1))*dt)...
          + start_dist(path(1));

%compute the intermediate gradients
dU(2:n-1) =  ((x_1-x_0).*(1-df_0*dt)-f_0*dt)/(noise_level^2*dt).*(1-df_0*dt)...
    -((x_2-x_1).*(1-df_1*dt)-f_1*dt)/(noise_level^2*dt).*(df2_1*dt.*(x_2-x_1)+1);

%compute the ending point gradient
dU(n) = ((path(n)-path(n-1))*(1-df(path(n-1))*dt)-f(path(n-1))*dt)/(noise_level^2*dt)*(1-df(path(n-1))*dt)...
          + end_dist(path(end));