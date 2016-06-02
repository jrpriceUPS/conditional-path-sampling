%clear all;close all;

addpath('../Functions')

k=10;

var_list = [1e-1,1e-2,1e-3,1e-4];
subplotdim = [2,2];
burn = 200;



SDE.initial = -1;

f{1}  = @(x) k*(-4*x.*(x.^2-1));
df{1} = @(x) k*(-12*x.^2+4);
    
drifts.f = f;
drifts.df = df;

cond.mean = 1;
cond.samples = 800;
cond.gap  = 0;
cond.burn = burn;

domain.dt       = 1e-3;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 0;

f = cell(1,1);
df = cell(1,1);

nLags = cond.samples-burn;




for i=1:length(noise_list)
    
    cond.var = var_list(i);
    
    
    output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);
    
    figure(2)
    subplot(subplotdim(1),subplotdim(2),i);
    correlates = autocorrelation(output.paths);
    plot(0:nLags-1,abs(correlates))
    title(sprintf('Variance = %g',cond.var))
    drawnow


    figure(3)
    subplot(subplotdim(1),subplotdim(2),i);
    plot(0:domain.dt:domain.endtime,output.paths)
    title(sprintf('Variacne = %g',cond.var))
    drawnow
 
    
end