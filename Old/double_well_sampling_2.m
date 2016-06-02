%clear all;close all;

addpath('../Functions')

%a double gaussian well potential with the particle starting in one of the
%wells

SDE.g1_height  =  5;                      %height of first gaussian
SDE.g1_mean    =  -1;                     %mean of first gaussian
SDE.g1_var     =  1/sqrt(2);              %variance of first gaussian

SDE.g2_height  =  20;                     %height of second gaussian
SDE.g2_mean    =  1;                      %mean of second gaussian
SDE.g2_var     =  0.25;                   %variance of second gaussian

SDE.left = -1;                            %left wall starts (for quartic boundary)
SDE.right = 1;                            %right wall starts (for quartic boundary)

SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = .01;
cond.samples = 10000;
cond.gap  = 50;
cond.burn = 50;

f  = @(x) SDE.g1_height*gaussian_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (4*(x-SDE.left).^3).*(x<SDE.left)...
    - (4*(x-SDE.right).^3).*(x>SDE.right);
df = @(x) SDE.g1_height*gaussian_2nd_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_2nd_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (12*(x-SDE.left).^2).*(x<SDE.left)...
    - (12*(x-SDE.right).^2).*(x>SDE.right);

drifts.f = cell(1,1);
drifts.f{1} = f;
drifts.df = cell(1,1);
drifts.df{1} = df;

domain.dt       = 1e-2;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 1;

output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);