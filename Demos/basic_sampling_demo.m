clear all; close all; clc
addpath('../Functions')

%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%

%degree of Brownian noise
SDE.noise = 1/2;



%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%

%depth of wells
k=10;

%standard quartic SDE drift and derivative of drift
f{1}  = @(x) k*(-4*x.*(x.^2-1));
df{1} = @(x) k*(-12*x.^2+4);

drifts.f = f;
drifts.df = df;



%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%

%time step
domain.dt       =  2^-7;   %timestep (resolution)
domain.endtime  =  1;      %end of simulation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional sampling parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mean and standard deviation of initial and final observation gaussians
cond.mean = [-1,1];  %start in left well, end in right
cond.std  = [.01,.01]; %tight distributions on each end

%functions relating to the potentials of initial and final observations
cond.start_neg_log    =  @(x) (x-cond.mean(1))^2/(2*cond.std(1)^2);
cond.end_neg_log      =  @(x) (x-cond.mean(2))^2/(2*cond.std(2)^2);
cond.start_d_neg_log  =  @(x) (x-cond.mean(1))/(cond.std(1)^2);
cond.end_d_neg_log    =  @(x) (x-cond.mean(2))/(cond.std(2)^2);

%draw an initial position from the initial distribution
cond.initial_pos      =  cond.std(1)*randn + cond.mean(1);

%number of samples
cond.samples = 1000;

%how infrequently to save paths (use when handling large numbers of
%possibly correlated data)
cond.gap  = 0;

%number of initial paths to disregard before recording results
cond.burn = 10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hybrid Monte Carlo parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time step of HMC

HMC_params.dt = 0;        %set to 0 for it to be automatically calculated
HMC_params.rnd_dt = 0.01; %amount of randomness in selecting HMC_params.dt (percentage)

%piecewise power law
if ~HMC_params.dt
    if domain.dt > 1/500
        HMC_params.dt = .2269*(domain.dt^.5434);
    else
        HMC_params.dt = .4632*(domain.dt^.6758);
    end
end

%number of time steps of HMC
HMC_params.L  = 0;    %set to 0 for it to be automatically calculated
HMC_params.rnd_L = 0; %amount of randomness in selecting HMC_params.L (percentage)

if ~HMC_params.L
    HMC_params.L  = round(5/HMC_params.dt);
end



%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%

%1 if plots should be generated, 0 otherwise
plots.show         =  1;   %1 if plots should be generated during simulation, 0 otherwise
plots.print_ratio  =  1;   %1 if acceptance rate should be printed after each step
plots.num_plotted  =  0;   %number of plots highlighted at end



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%Actual simulation%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%


output = conditional_path(SDE,drifts,cond,domain,HMC_params,plots);
output.accept_rate



%%%%%%%%%%%%%%%%%
%Post processing%
%%%%%%%%%%%%%%%%%

%compute the number of lag times we can compute with our Fourier method
nLags = (cond.samples-cond.burn)/2;


paths_direct=figure(1);
plot(0:domain.dt:domain.endtime,output.paths)
hold on

if plots.num_plotted
    plot(0:domain.dt:domain.endtime,output.paths(:,end:-1*ceil((cond.samples - cond.burn)/plots.num_plotted):1),'k-','LineWidth',2)
end

title('Sampled Paths','FontSize',16)
axis([0,1,-2,2])

trans_times_direct=figure(2);
[m,n]=size(output.paths);
plot(1:n,sum(output.paths<0)/m)
title('Transition Times','FontSize',16)
axis([1,n,0,1])

autocorr=figure(3);
correlates = autocorrelation(output.paths);
plot(0:nLags-1,correlates(:,1:nLags))
title('Autocorrelation Function (all times)','FontSize',16)
axis([0,nLags-1,-1,1])

%save plots
%saveas(paths_direct,'paths_direct.png')
%saveas(trans_times_direct,'trans_times_direct.png')
%saveas(autocorr,'autocorrelation_direct_all.png')

%save animation to gif
% figure,
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% 
% for k=1:cond.samples - cond.burn,
%     plot(0:domain.dt:domain.endtime,output.paths(:,1:k),'-');
%     f=getframe;
%     im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
% end
% 
% imwrite(im,map,'pathsAnimation.gif','DelayTime',0,'LoopCount',inf);
