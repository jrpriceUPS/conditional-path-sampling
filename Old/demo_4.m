clear all;close all;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metadynamical parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%

meta.weight    =  1;         %initial weight of gaussian deposits
meta.width     =  0.2;       %width of deposits
meta.tempered  =  1;         %logical variable (1 if well-tempered, 0 otherwise)
meta.freq      =  50;        %how many timesteps to take between deposits
meta.temp      =  100/kB;    %temperature of system (for well-tempered)

domain_meta.dt        =  0.001;  %timestep
domain_meta.endtime   =  1;      %end time
domain_meta.periodic  =  0;      %logical variable (1 if periodic domain, 0 otherwise)


plots_meta.show        =  0;             %logical variable (1 if metadynamics should be plotted, 0 otherwise)
plots_meta.axes        =  [-4,4,-20,4];  %axes of metadynamics
plots_meta.resolution  =  1000;          %x-axis resolution of metadynamics



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional path dynamics parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cond.mean = 1;        %mean of end state on which we are conditioning      
cond.var = .01;       %variance of end state on which we are conditioning
cond.samples = 100;   %number of samples at each level of detail
cond.n_levels = 10;   %number of levels of detail

domain_cond.dt       = 0.001;                 %timestep
domain_cond.endtime  = 1;                     %end time
domain_cond.periodic = domain_meta.periodic;  %logical variable (1 if periodic domain, 0 otherwise)

plots_cond.show = 1;    %logical variable (1 if paths displayed after each level, 0 otherwise)
plots_cond.gap  = 20;   %gap between paths (empirically determine gap needed for independence)
plots_cond.burn = 5;    %number of paths to ignore at start of simulation

HMC_params.dt = 0.005;  %timestep of hybrid monte carlo step
HMC_params.T  = 1;      %end time of each hybrid monte carlo step



%%%%%%%%%%%%%%%%%%%%
%Actual simulations%
%%%%%%%%%%%%%%%%%%%%

data = metadynamics(SDE,meta,domain_meta,plots_meta);
drifts = meta_processing(cond,SDE,data,meta);
paths = conditional_path_v2(SDE,drifts,cond,domain_cond,HMC_params,plots_cond);