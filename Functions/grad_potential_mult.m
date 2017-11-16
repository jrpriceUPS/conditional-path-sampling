function dU = grad_potential_mult(path,start_dist,end_dist,f,df,noise_level,dt)
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

%compute the resolution of the path
M = size(path,1);

%compute the dimension
n = size(path,2);

%initialize result
dU = zeros(M,n);

x = zeros(M,n); %x(i,j) = f(path(i,:))[j] i.e. the jth component in f(x_i) where x_i is the ith point on the path (the ith row in path)
for i = 1:M
    x(i,:) = f(path(i,:));
end

H = zeros(M,n); %rows are the diagonals of the Hessian of the Mth point on the path
for i = 1:M %for each point on the path
    H(i,:) = diag(df(path(i,:))).'; %ith row is the diagonals of the Hessian
end

for j = 1:n
    
    %start point
    dU(1,j) = -(path(2,j)-(path(1,j)+x(1,j)*dt))/(noise_level^2*dt).*(1+H(1,j)*dt);
    
    %intermediate points
    dU(2:M-1,j) = (path(2:M-1,j)-(path(1:M-2,j)+x(1:M-2,j)*dt))/(noise_level^2*dt)...
            - (path(3:M,j)-(path(2:M-1,j)+x(2:M-1,j)*dt))/(noise_level^2*dt).*(1+H(2:M-1,j)*dt);
    
    %end point
    dU(M,j) = (path(M) - (path(M-1,j)+x(M-1,j)*dt))/(noise_level^2*dt);

end

dU(1,:) = dU(1,:) + start_dist(path(1,:)).';
dU(end,:) = dU(end,:) + end_dist(path(end,:)).';