clear all;close all;

addpath('../Functions')

k=10;

SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = .1;
cond.samples = 200;

f = cell(1,1);
df = cell(1,1);

f{1}  = @(x) k*(-4*x.*(x.^2-1));
df{1} = @(x) k*(-12*x.^2+4);

drifts.f = f;
drifts.df = df;

domain.dt       = 1e-3;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 0;

output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);



SDE.initial_path  =  output.paths(:,end);
cond.samples = 100;
cond.var = 0.5;
plots.show=1;

output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);


