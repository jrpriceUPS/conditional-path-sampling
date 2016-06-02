clear all;close all;

addpath('../Functions')
k=10;

k_list = [2.5,5,7.5,10];
subplotdim = [4,1];
burn = 200;


SDE.noise = 1/2;
SDE.initial = -1;

cond.mean = 1;
cond.var = .01;
cond.samples = 100000;
cond.gap = 0;
cond.burn = burn;

domain.dt       = 0.001;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 0;

f = cell(1,1);
df = cell(1,1);

nLags = (cond.samples-burn)/2;


corr = figure(1);
paths = figure(2);
trans_times = figure(3);


for i=1:length(k_list)
    
    k = k_list(i);
    
    f{1}  = @(x) k*(-4*x.*(x.^2-1));
    df{1} = @(x) k*(-12*x.^2+4);
    
    drifts.f = f;
    drifts.df = df;
    
    
    output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);
    
    figure(1);
    subplot(subplotdim(1),subplotdim(2),i);
    correlates = autocorrelation(output.paths(2:end,:));
    plot(0:nLags-1,correlates)
    title(sprintf('k = %g',k))
    
    figure(2);
    subplot(subplotdim(1),subplotdim(2),i);
    plot(0:domain.dt:domain.endtime,output.paths)
    title(sprintf('k = %g',k))
    
    figure(3);
    subplot(subplotdim(1),subplotdim(2),i);
    plot(1:cond.samples-burn,sum(output.paths<0)/1001);
    axis([1,cond.samples-burn,0,1])
    title(sprintf('k = %g',k))
 
    
end

saveas(corr,'corr.png')
saveas(paths,'paths.png')
saveas(trans_times,'trans_times.png')