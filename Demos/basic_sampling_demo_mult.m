clear all; close all; clc
addpath('../Functions')


%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%

%degree of Brownian noise
SDE.noise = 1/2;


%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%

%depth of wells
k=5;

%standard quartic SDE drift and derivative of drift
%f = -grad(U), where U is the potential energy landscape
%here the input x is a row vector of length 2
%U     = @(x) -2*k*x(1).^2 + x(2).^2 + (x(1).^2).*(x(2).^2) + k*(x(1).^4) + x(2).^4;
f{1}  = @(x) [4*k*x(1) - 2*x(1).*(x(2).^2) - 4*k*(x(1).^3), -2*x(2) - 2*(x(1).^2).*x(2) - 4*(x(2).^3)];
df{1} = @(x) [4*k - 2*(x(2).^2) - 12*k*(x(1).^2), -4*(x(1).*x(2));
              -4*(x(1).*x(2)), -2 - 2*(x(1).^2) - 12*(x(2).^2)]; % negative Hessian of U

drifts.f = f;
drifts.df = df;


%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%

%time step
domain.dt       =  1/300;  %timestep (resolution)
domain.endtime  =  1;      %end of simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional sampling parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mean and standard deviation of initial and final observation gaussians
cond.mean = [-1 0; 1 0];  %first column is x coordiante, second column is y coordiante. Top row is start point, bottom row is end point.
cond.std  = [.01,.01];    %tight distributions on each end

n = 2; %dimension of the paths
cond.initial_pos      =  zeros(1,n);
for i = 1:n
    cond.initial_pos(i) = cond.std(1)*randn + cond.mean(1,i);
end

sigma  =  SDE.noise;

dt = domain.dt;
T  = domain.endtime;

X0  =  cond.initial_pos;

%initialize the time domain
t = 0:dt:T;

%compute initial paths with Langevin dynamics
initial_path       =  zeros(length(t),n); %first column is x coordinates, second column is y coordinates; each row is a point in R^2
initial_path(1,:)  =  X0;
    
for i = 1:(T/dt)
        
    initial_path(i+1,:) = initial_path(i,:) + f{1}(initial_path(i,:))*dt + sigma*sqrt(dt)*randn(1,n);
    
end

%prints the initial path
%initial_path

%plotting the free energy landscape (heatmap)
[X,Y] = meshgrid(-2:.1:2);                                
Z = @(a) -2*a*(X.^2) + (Y.^2) + (X.^2).*(Y.^2) + a*(X.^4) + Y.^4;
surf(X,Y,Z(k))
view(2)

hold on

%plotting the initial path
plot(initial_path(:,1)', initial_path(:,2)', 'r')
axis([-2,2,-2,2])