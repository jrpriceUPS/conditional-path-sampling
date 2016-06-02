clear all;close all;

addpath('../Functions')

SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = [-1,1];
cond.std = [.1,.1];

cond.start_neg_log    =  @(x) (x-cond.mean(1))^2/(2*cond.std(1)^2);
cond.end_neg_log      =  @(x) (x-cond.mean(2))^2/(2*cond.std(2)^2);
cond.start_d_neg_log  =  @(x) (x-cond.mean(1))/(cond.std(1)^2);
cond.end_d_neg_log    =  @(x) (x-cond.mean(2))/(cond.std(2)^2);
cond.initial_pos      =  cond.std(1)*randn + cond.mean(1);

cond.samples = 1000;
cond.gap  = 0;
cond.burn = 20;

drift_mod_list = [7.5,8,8.5,9,9.5,10];

f = cell(length(drift_mod_list),1);
df = cell(length(drift_mod_list),1);

f_0  = @(x) -4*x.*(x.^2-1);
df_0 = @(x) -12*x.^2+4;

for i=1:length(drift_mod_list)
    f{i} = @(x)  drift_mod_list(i)*f_0(x);
    df{i} = @(x) drift_mod_list(i)*df_0(x);
end

drifts.f = f;
drifts.df = df;

domain.dt       = 1e-3;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 1;
plots.subplot_dim = [1,1];

output = conditional_path(SDE,drifts,cond,domain,HMC_params,plots);


nLags = (cond.samples-cond.burn)/2;

corr=figure(1);
paths=figure(2);
skip=50;


figure(1)

correlates = autocorrelation(output.paths);
plot(0:nLags-1,correlates)
title(sprintf('k = %g',drift_mod_list(i)*10))
axis([0,nLags-1,-1,1])

figure(2)
plot(0:domain.dt:domain.endtime,output.paths)
title(sprintf('k = %g',drift_mod_list(i)*10))
axis([0,1,-2,2])
