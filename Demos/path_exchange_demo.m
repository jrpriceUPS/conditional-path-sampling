clear all;close all;
addpath('../Functions')

%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%

%degree of Brownian noise
SDE.noise = 1/2;



%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%

%list of multipliers of standard quartic SDE potential (higher = deeper
%wells)
drift_mod_list = [7.5,8,8.5,9,9.5,10];

%initialize drift function cells
f = cell(length(drift_mod_list),1);
df = cell(length(drift_mod_list),1);

%standard quartic SDE drift and derivative of drift
f_0  = @(x) -4*x.*(x.^2-1);
df_0 = @(x) -12*x.^2+4;

%save multiplied drifts
for i=1:length(drift_mod_list)
   f{i} = @(x)  drift_mod_list(i)*f_0(x);
   df{i} = @(x) drift_mod_list(i)*df_0(x);
end

drifts.f = f;
drifts.df = df;



%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%

domain.dt       =  2^-10;  %timestep
domain.endtime  =  1;      %end of simulation



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
cond.initial_pos      =  cond.std(1)*randn(1,length(drift_mod_list)) + cond.mean(1);

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
HMC_params.dt = 0.005;

%end time of HMC
HMC_params.T  = 1;



%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%

plots.show = 0;            %1 if plots should be generated during simulation, 0 otherwise
plots.subplot_dim = [3,2]; %dimensions of subplot array
plots.num_plotted  =  10;  %number of plots highlighted at end



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%Actual simulation%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%


output = conditional_path_parallel(SDE,drifts,cond,domain,HMC_params,plots);



%%%%%%%%%%%%%%%%%
%Post processing%
%%%%%%%%%%%%%%%%%

%compute the number of lag times we can compute with our Fourier method
nLags = (cond.samples-cond.burn)/2;

%plot paths, highlight a fraction of them
paths_direct=figure(1);
plot(0:domain.dt:domain.endtime,squeeze(output.paths(:,end,:)))
hold on
plot(0:domain.dt:domain.endtime,squeeze(output.paths(:,end,cond.samples/plots.num_plotted:cond.samples/plots.num_plotted:end)),'k-','LineWidth',2)
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