function U = potential_mult(path,start_dist,end_dist,f,noise_level,dt)
%Takes the path and computes the potential U under the explicit
%Euler-Maruyama scheme
%
%INPUTS:
%
%path         =  the path; an Mxn matrix; each row is a point
%start_dist   =  the negative logarithm of the starting distribution
%end_dist     =  the negative logarithm of the ending distribution
%f            =  the drift function
%noise_level  =  the level of Brownian noise (sigma)
%dt           =  the time increment

%compute the potential (see May note)

x = zeros(size(path,1),size(path,2)); %x(i,j) = f(path(i,:))[j] i.e. the jth component in f(x_i) where x_i is the ith point on the path (the ith row in path)
for i = 1:size(path,1)
    x(i,:) = f(path(i,:));
end

U = 0;
for j = 1:size(path,2)
    U = U + sum((path(2:end,j)-(path(1:end-1,j)+x(1:end-1,j)*dt)).^2/(2*noise_level^2*dt));
end

U = U + start_dist(path(1,:)) + end_dist(path(end,:));