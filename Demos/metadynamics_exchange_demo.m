clear all;close all;

addpath('../Functions')
addpath('../Metadynamics')

%physical constant
kB = 1.38064852*10^-23;


%%%%%%%%%
%The SDE%
%%%%%%%%%

SDE.depth = 10;    %depth of wells

SDE.drift        =  @(x) SDE.depth*-4*x.*(x.^2-1);  %drift
SDE.drift_deriv  =  @(x) SDE.depth*(-12*x.^2+4);    %derivative of drift
SDE.potential    =  @(x) SDE.depth*(x.^4-2*x.^2);   %potential
SDE.noise        =  1/2;                            %noise level
SDE.initial      =  -1;                             %initial condition



%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metadynamical parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%

meta.weight    =  0.5;       %initial weight of gaussian deposits
meta.width     =  0.5;       %width of deposits
meta.tempered  =  1;         %logical variable (1 if well-tempered, 0 otherwise)
meta.freq      =  10;        %how many timesteps to take between deposits
meta.temp      =  1000/kB;   %temperature of system (for well-tempered)

domain_meta.dt        =  0.01;  %timestep
domain_meta.endtime   =  2;     %end time

plots_meta.show        =  0;                                %logical variable (1 if metadynamics should be plotted, 0 otherwise)
plots_meta.axes        =  [-2,2,-SDE.depth,SDE.depth*3/2];  %axes of metadynamics
plots_meta.resolution  =  1000;                             %x-axis resolution of metadynamics



%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%

domain_cond.dt       =  2^-10;
domain_cond.endtime  =  1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional sampling parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of levels of intermediate potentials (excluding base potential)
cond.n_levels = 3;

%mean and standard deviation of initial and final observation gaussians
cond.mean = [-1,1];  %start in left well, end in right
cond.std  = [.1,.1]; %tight distributions on each end

%functions relating to the potentials of initial and final observations
cond.start_neg_log    =  @(x) (x-cond.mean(1))^2/(2*cond.std(1)^2);
cond.end_neg_log      =  @(x) (x-cond.mean(2))^2/(2*cond.std(2)^2);
cond.start_d_neg_log  =  @(x) (x-cond.mean(1))/(cond.std(1)^2);
cond.end_d_neg_log    =  @(x) (x-cond.mean(2))/(cond.std(2)^2);

%draw an initial position from the initial distribution
cond.initial_pos      =  cond.std(1)*randn(1,cond.n_levels+1) + cond.mean(1);

%number of samples
cond.samples = 1000;

%how infrequently to save paths (use when handling large numbers of
%possibly correlated data)
cond.gap  = 0;

%number of initial paths to disregard before recording results
cond.burn = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hybrid Monte Carlo parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time step of HMC
HMC_params.dt = 0.005;

%end time of HMC
HMC_params.T  = 1;



%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%


plots_cond.show = 0;  %1 if plots should be generated during simulation, 0 otherwise
plots_cond.subplot_dim = [2,2];  %dimensions of subplot array
plots_cond.num_plotted  =  10;  %number of plots highlighted at end


%%%%%%%%%%%%%%%%%%%%
%Actual simulations%
%%%%%%%%%%%%%%%%%%%%

data = metadynamics(SDE,meta,domain_meta,plots_meta);
drifts = meta_processing(cond,SDE,data,meta);
output = conditional_path_parallel(SDE,drifts,cond,domain_cond,HMC_params,plots_cond);



%%%%%%%%%%%%%%%%%
%Post processing%
%%%%%%%%%%%%%%%%%

%compute the number of lag times we can compute with our Fourier method
nLags = (cond.samples-cond.burn)/2;

%plot paths, highlight a fraction of them
paths_direct=figure(1);
plot(0:domain_cond.dt:domain_cond.endtime,squeeze(output.paths(:,end,:)))
hold on
plot(0:domain_cond.dt:domain_cond.endtime,squeeze(output.paths(:,end,cond.samples/plots_cond.num_plotted:cond.samples/plots_cond.num_plotted:end)),'k-','LineWidth',2)
title('Sampled Paths','FontSize',16)
axis([0,1,-2,2])

%heuristically plot the time of transition for each path
trans_times_direct=figure(2);
[m,n]=size(squeeze(output.paths(:,end,:)));
plot(1:n,sum(squeeze(output.paths(:,end,:))<0)/m)
title('Transition Times','FontSize',16)
axis([1,n,0,1])

%compute and plot autocorrelations
autocorr=figure(3);
correlates = autocorrelation(squeeze(output.paths(:,end,:)));
plot(0:nLags-1,correlates(:,1:nLags))
title('Autocorrelation Function (all times)','FontSize',16)
axis([0,nLags-1,-1,1])

%save plots
saveas(paths_direct,'paths_direct.png')
saveas(trans_times_direct,'trans_times_direct.png')
saveas(autocorr,'autocorr.png')