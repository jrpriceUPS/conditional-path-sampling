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

dt_max_res         =  -10;  %two to this power is highest level of resolution
dt_num_lev         =  4;    %number of resolutions (goes up in powers of two)
domain.endtime     =  1;    %end of simulation

%(row) vector of resolutions
domain.dt = 2.^(dt_num_lev-(1:dt_num_lev)+dt_max_res);



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
cond.initial_pos      =  cond.std(1)*randn(1,dt_num_lev) + cond.mean(1);

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
HMC_params.dt = []; %set to an empty or zero vector for it to be automatically calculated

if (length(HMC_params.dt) ~= dt_num_lev) || (max(HMC_params.dt) == 0)
    HMC_params.dt = zeros(1, dt_num_lev);
    for j = 1:dt_num_lev 
        if domain.dt(j) > 1/500
            HMC_params.dt(j) = .2269*(domain.dt(j)^.5434);
        else
            HMC_params.dt(j) = .4632*(domain.dt(j)^.6758);
        end
    end
end

%number of HMC steps
HMC_params.L = round(4.75./HMC_params.dt);



%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%

plots.show         = 1;      %1 if plots should be generated during simulation, 0 otherwise
plots.print_ratio  = 1;      %1 if acceptance rate should be printed after each step, 0 otherwise
plots.num_plotted  = 10;     %number of plots highlighted at end
plots.subplot_dim  = [ceil(dt_num_lev / 2),2];  %dimensions of subplot array



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
plot(0:2^dt_max_res:domain.endtime,output.paths{dt_num_lev})
hold on

if plots.num_plotted
    plot(0:domain.dt(end):domain.endtime,output.paths{dt_num_lev}(:,end:-1*ceil((cond.samples - cond.burn)/plots.num_plotted):1),'k-','LineWidth',2)
end

title('Sampled Paths','FontSize',16)
axis([0,1,-2,2])

%heuristically plot the time of transition for each path
trans_times_direct=figure(2);
[m,n]=size(output.paths{dt_num_lev});
plot(1:n,sum(output.paths{dt_num_lev}<0)/m)
title('Transition Times','FontSize',16)
axis([1,n,0,1])

%compute and plot autocorrelations
autocorr=figure(3);
correlates = autocorrelation(output.paths{dt_num_lev});
plot(0:nLags-1,correlates(:,1:nLags))
title('Autocorrelation Function (all times)','FontSize',16)
axis([0,nLags-1,-1,1])

%save plots
% saveas(paths_direct,'paths_direct.png')
% saveas(trans_times_direct,'trans_times_direct.png')
% saveas(autocorr,'autocorr.png')