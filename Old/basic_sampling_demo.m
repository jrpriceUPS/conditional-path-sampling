%clear all;close all;

addpath('../Functions')

SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = .01;
cond.samples = 100000;
cond.gap  = 0;
cond.burn = 0;

f  = @(x) 10*(-4*x.*(x.^2-1));
df = @(x) 10*(-12*x.^2+4);

drifts.f = cell(1,1);
drifts.f{1} = f;
drifts.df = cell(1,1);
drifts.df{1} = df;

domain.dt       = 1e-3;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 0;

output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);