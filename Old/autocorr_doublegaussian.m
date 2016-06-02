clear all;close all;

addpath('../Functions')

burn = 5;


SDE.noise = 1/2;
SDE.initial = -1;

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



cond.mean = 1;
cond.var = .01;
cond.samples = 500;

domain.dt       = 1e-3;
domain.endtime  = 1;
domain.periodic = 0;

HMC_params.dt = 0.005;
HMC_params.T  = 1;

plots.show = 1;
plots.gap  = 150;
plots.burn = 150;

f = cell(1,1);
df = cell(1,1);

nLags = cond.samples-burn;

correlates = zeros(domain.endtime/domain.dt+1,nLags+1);





f{1}  = @(x) SDE.g1_height*gaussian_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (4*(x-SDE.left).^3).*(x<SDE.left)...
    - (4*(x-SDE.right).^3).*(x>SDE.right);
df{1} = @(x) SDE.g1_height*gaussian_2nd_deriv(x,SDE.g1_mean,SDE.g1_var,1)...
    + SDE.g2_height*gaussian_2nd_deriv(x,SDE.g2_mean,SDE.g2_var,1)...
    - (12*(x-SDE.left).^2).*(x<SDE.left)...
    - (12*(x-SDE.right).^2).*(x>SDE.right);

drifts.f = f;
drifts.df = df;


output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots);

output.paths = output.paths(:,burn+1:end);

figure(2)
correlates(:,:) = autocorrelation(nLags,output.paths);
plot(0:nLags,correlates)
axis([0,nLags,-1,1])
title('Autocorrelation')
drawnow
