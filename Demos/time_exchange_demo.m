clear all;close all;
addpath('../Functions')

%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%

%degree of Brownian noise
SDE.noise = 1;



%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%

%depth of wells
k=10;

%standard quartic SDE drift and derivatives of drift
f{1}   = @(x) k*(-4*x.*(x.^2-1));
df{1}  = @(x) k*(-12*x.^2+4);
df2{1} = @(x) k*(-24*x);

drifts.f   = f;
drifts.df  = df;
drifts.df2 = df2;



%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%

domain.dt_max_res  =  -10;  %two to this power is highest level of resolution
domain.dt_num_lev  =  4;    %number of resolutions (go up in powers of two)
domain.endtime     =  1;    %end of simulation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional sampling parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mean and standard deviation of initial and final observation gaussians
cond.mean = [-1,1];  %start in left well, end in right
cond.std  = [.1,.1]; %tight distributions on each end

%functions relating to the potentials of initial and final observations
cond.start_neg_log    =  @(x) (x-cond.mean(1))^2/(2*cond.std(1)^2);
cond.end_neg_log      =  @(x) (x-cond.mean(2))^2/(2*cond.std(2)^2);
cond.start_d_neg_log  =  @(x) (x-cond.mean(1))/(cond.std(1)^2);
cond.end_d_neg_log    =  @(x) (x-cond.mean(2))/(cond.std(2)^2);

%draw an initial position from the initial distribution
cond.initial_pos      =  cond.std(1)*randn(1,domain.dt_num_lev) + cond.mean(1);

%number of samples
cond.samples = 1000;

%how infrequently to save paths (use when handling large numbers of
%possibly correlated data)
cond.gap  = 0;

%number of initial paths to disregard before recording results
cond.burn = 10;

%determine whether to use explicit (=0) or implicit Euler-Maruyama (=1)
cond.implicit = 0;



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

plots.show         = 0;      %1 if plots should be generated during simulation, 0 otherwise
plots.subplot_dim  = [2,2];  %dimensions of subplot array
plots.num_plotted  =  10;    %number of plots highlighted at end



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%Actual simulation%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%


output = conditional_path_time(SDE,drifts,cond,domain,HMC_params,plots);



%%%%%%%%%%%%%%%%%
%Post processing%
%%%%%%%%%%%%%%%%%

%compute the number of lag times we can compute with our Fourier method
nLags = (cond.samples-cond.burn)/2;


%plot paths, highlight a fraction of them
paths_direct=figure(1);
plot(0:2^domain.dt_max_res:domain.endtime,output.paths{4})
hold on
plot(0:2^domain.dt_max_res:domain.endtime,output.paths{4}(:,cond.samples/plots.num_plotted:cond.samples/plots.num_plotted:end),'k-','LineWidth',2)
title('Sampled Paths','FontSize',16)
axis([0,1,-2,2])

%heuristically plot the time of transition for each path
trans_times_direct=figure(2);
[m,n]=size(output.paths{4});
plot(1:n,sum(output.paths{4}<0)/m)
title('Transition Times','FontSize',16)
axis([1,n,0,1])

%compute and plot autocorrelations
autocorr=figure(3);
correlates = autocorrelation(output.paths{4});
plot(0:nLags-1,correlates(:,1:nLags))
title('Autocorrelation Function (all times)','FontSize',16)
axis([0,nLags-1,-1,1])

%save plots
saveas(paths_direct,'paths_direct.png')
saveas(trans_times_direct,'trans_times_direct.png')
saveas(autocorr,'autocorr.png')