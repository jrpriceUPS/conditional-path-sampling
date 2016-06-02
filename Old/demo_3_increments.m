%clear all;close all;

addpath('../Functions')

n_levels=10;

SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = 0.01;
cond.samples = 100;

d_drift = 1/n_levels;

drift_mod_list = d_drift:d_drift:1;

f = cell(length(drift_mod_list),1);
df = cell(length(drift_mod_list),1);

f_0  = @(x) 10*(-4*x.*(x.^2-1));
df_0 = @(x) 10*(-12*x.^2+4);

for i=1:length(drift_mod_list)
   f{i} = @(x)  drift_mod_list(i)*f_0(x);
   df{i} = @(x) drift_mod_list(i)*df_0(x);
   %f{i} = f_0;
   %df{i} = df_0;
end

drifts.f = f;
drifts.df = df;

domain.dt       = 0.01;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.006275;
HMC_params.T  = 1;

plots.show = 1;

output = conditional_path(SDE,drifts,cond,domain,HMC_params,plots);