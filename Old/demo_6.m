%clear all;close all;

addpath('../Functions')

%physical constant
kB = 1.38064852*10^-23;


%%%%%%%%%
%The SDE%
%%%%%%%%%

%a double gaussian well potential with the particle starting in one of the
%wells

SDE.g1_height  =  5;                      %height of first gaussian
SDE.g1_mean    =  1;                      %mean of first gaussian
SDE.g1_var     =  1/sqrt(2);              %variance of first gaussian

SDE.g2_height  =  20;                     %height of second gaussian
SDE.g2_mean    =  -1;                     %mean of second gaussian
SDE.g2_var     =  0.25;                   %variance of second gaussian

SDE.left = -1;                            %left wall starts (for quartic boundary)
SDE.right = 1;                            %right wall starts (for quartic boundary)

%drift, derivative of drift, and potential for double well with quartic
%boundaries
SDE.drift        =  @(x) SDE.g1_height*gaussian_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (4*(x-SDE.left).^3).*(x<SDE.left)...
    - (4*(x-SDE.right).^3).*(x>SDE.right);
SDE.drift_deriv  =  @(x) SDE.g1_height*gaussian_2nd_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_2nd_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (12*(x-SDE.left).^2).*(x<SDE.left)...
    - (12*(x-SDE.right).^2).*(x>SDE.right);
SDE.potential    =  @(x) -SDE.g1_height*gaussian(x,SDE.g1_mean,SDE.g1_var,1)...
    - SDE.g2_height*gaussian(x,SDE.g2_mean,SDE.g2_var,1)...
    + ((x-SDE.left).^4).*(x<SDE.left)...
    + ((x-SDE.right).^4).*(x>SDE.right);

SDE.noise      =  1;              %noise level
SDE.initial    =  SDE.g2_mean;    %initial condition



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional path dynamics parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cond.mean = 1;        %mean of end state on which we are conditioning
cond.var = .01;       %variance of end state on which we are conditioning
cond.samples = 1000;   %number of samples at each level of detail
cond.n_levels = 0;    %number of levels of detail

domain_cond.dt       = 0.001;                 %timestep
domain_cond.endtime  = 1;                     %end time
domain_cond.periodic = 0;  %logical variable (1 if periodic domain, 0 otherwise)

plots_cond.show = 0;    %logical variable (1 if paths displayed after each level, 0 otherwise)

HMC_params.dt = 0.005;  %timestep of hybrid monte carlo step
HMC_params.T  = 1;      %end time of each hybrid monte carlo step

num_runs = 360;
path_list = zeros(domain_cond.endtime/domain_cond.dt+1,num_runs);
accept_list = zeros(num_runs,1);
t = 0:domain_cond.dt:domain_cond.endtime;

%%%%%%%%%%%%%%%%%%%%
%Actual simulations%
%%%%%%%%%%%%%%%%%%%%

drifts.f   =  cell(1,1);
drifts.f{1}   =  SDE.drift;
drifts.df  =  cell(1,1);
drifts.df{1}  =  SDE.drift_deriv;

for j=1:num_runs
    paths = conditional_path_v2(SDE,drifts,cond,domain_cond,HMC_params,plots_cond);
    accept_list(j) = mean(paths.accept_rate);
    path_list(:,j) = paths.paths(:,end);
end

plot(t,path_list);