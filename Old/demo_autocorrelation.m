%clear all;close all;

addpath('../Functions')
k=10;

dt_list = [0.5,0.1,0.01,0.001];
subplotdim = [3,1];
burn = 200;


SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = .01;
cond.samples = 1000;
cond.gap = 0;
cond.burn = burn;

domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 0;

f = cell(1,1);
df = cell(1,1);
df2 = cell(1,1);

nLags = (cond.samples-burn)/2;





for i=1:length(dt_list)
    
    domain.dt = dt_list(i);
    
    f{1}   = @(x) k*(-4*x.*(x.^2-1));
    df{1}  = @(x) k*(-12*x.^2+4);
    df2{1} = @(x) k*(-24*x);
    
    drifts.f = f;
    drifts.df = df;
    drifts.df2 = df2;
    
    
    output = conditional_path(SDE,drifts,cond,domain,HMC_params,plots);
    
    figure(1)
    subplot(subplotdim(1),subplotdim(2),i);
    correlates = autocorrelation(output.paths(2:end,:));
    plot(0:nLags-1,correlates)
    title(sprintf('dt = %g',domain.dt))
    drawnow
    
    figure(2)
    subplot(subplotdim(1),subplotdim(2),i);
    plot(0:domain.dt:domain.endtime,output.paths)
    title(sprintf('dt = %g',domain.dt))
    drawnow
 
    
end