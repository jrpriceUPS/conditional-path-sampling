clear all; close all; clc
addpath('../Functions')
tic
%%%%%%%%%%%
%Dimension%
%%%%%%%%%%%
n = 2;
%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%
%degree of Brownian noise
SDE.noise = 1;
%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%
%depth of wells
k = 1;
%standard quartic SDE drift and derivative of drift
%f = -grad(U), where U is the potential energy landscape
%the input x is a row vector of length n
U     = @(x) -2*k*x(1)^2 + x(2)^2 + (x(1)^2)*(x(2)^2) + k*(x(1)^4) + x(2)^4;
f  = @(x) [4*k*x(1) - 2*x(1)*(x(2)^2) - 4*k*(x(1)^3), -2*x(2) - 2*(x(1)^2)*x(2) - 4*(x(2)^3)];
df = @(x) [4*k - 2*(x(2)^2) - 12*k*(x(1)^2), -4*x(1)*x(2);
              -4*x(1)*x(2), -2 - 2*(x(1)^2) - 12*(x(2)^2)]; % negative Hessian of U

drifts.f = f;
drifts.df = df;
drifts.U = U;
drifts.k = k;

%%%%%%%%%%%%%%%%%%%
%Domain parameters%
%%%%%%%%%%%%%%%%%%%
dt_max_res         =  -7;  %two to this power is highest level of resolution
dt_num_lev         =  4;    %number of resolutions (goes up in powers of two)
domain.endtime     =  1;    %end of simulation

%(row) vector of resolutions
domain.dt = 2.^(dt_num_lev-(1:dt_num_lev)+dt_max_res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conditional sampling parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mean and standard deviation of initial and final observation gaussians
cond.mean       = [-1 0; 1 0];        %top row is start point, bottom row is end point.
cond.cov_start  = [.0001 0; 0 .0001];   %(invertible) covariance matrices of start
cond.cov_end    = [.0001 0; 0 .0001];   %and end distributions
cond.std        = [.01 .01];
%0 correlation, standard deviation of .01 in each variable

%functions relating to the potentials of initial and final observations
%negative log of a multivariate normal distribution
%this part needs to be optimized
cond.start_neg_log    =  @(x) (x - cond.mean(1,:))*inv(cond.cov_start)*(x - cond.mean(1,:)).'/2;
cond.end_neg_log      =  @(x) (x - cond.mean(2,:))*inv(cond.cov_end)*(x - cond.mean(2,:)).'/2;
cond.start_d_neg_log  =  @(x) (inv(cond.cov_start) + inv(cond.cov_start)')*(x - cond.mean(1,:)).'/2;
cond.end_d_neg_log    =  @(x) (inv(cond.cov_end) + inv(cond.cov_end)')*(x - cond.mean(2,:)).'/2;

%draw an initial position from the initial distribution
cond.initial_pos = zeros(1,2,dt_num_lev);
for i = 1:dt_num_lev
cond.initial_pos(:,:,i) = cond.mean(1,:) + cond.std(1)*randn(1,n);
end

%the exchange frequency, which means the program exchage a pair of path 
%every cond.exchange_frequency samples  
cond.exchange_frequency = 1;

%number of samples
cond.samples = 2000;

%how infrequently to save paths (use when handling large numbers of
%possibly correlated data)
cond.gap  = 0;

%number of initial paths to disregard before recording results
cond.burn = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hybrid Monte Carlo parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time step of HMC
HMC_params.dt = [];        %set to 0 for it to be automatically calculated
HMC_params.rnd_dt = 0.01; %amount of randomness in selecting HMC_params.dt (percentage)

%piecewise power law
if  (length(HMC_params.dt) ~= dt_num_lev) || (max(HMC_params.dt) == 0)
    for j = 1:dt_num_lev 
        if domain.dt(j) > 1/500
            HMC_params.dt(j) = .2269*(domain.dt(j)^.5434);
        else
            HMC_params.dt(j) = .4632*(domain.dt(j)^.6758);
        end
    end
end

%number of time steps of HMC
HMC_params.rnd_L = 0; %amount of randomness in selecting HMC_params.L (percentage)
HMC_params.L  = round(5./HMC_params.dt); %automatically calculated from the dt

%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%

plots.show         = 0;      %1 if plots should be generated during simulation, 0 otherwise
plots.print_ratio  = 0;      %1 if acceptance rate should be printed after each step, 0 otherwise
plots.num_plotted  = 0;     %number of plots highlighted at end
plots.subplot_dim  = [ceil(dt_num_lev / 2),2];  %dimensions of subplot array

%%%%%%%%%%%%%%%%%%%%%
%%Actual simulation%%
%%%%%%%%%%%%%%%%%%%%%

output = conditional_path_time_mult(SDE,drifts,cond,domain,HMC_params,plots,n);
%output.accept_rate{i} is the accept rate of the ith resolution
%output.path{i} is the path of the ith resolution 

%%%%%%%%%%%%%%%%%
%Post processing%
%%%%%%%%%%%%%%%%%

%make a plot of transitioon time
%initialize a figure to hold the transition time of the highest resolution
trans_times_direct=figure(2);
%save the number of points in M
M = size(output.paths{4},1);
%save the number of paths in num_paths
num_paths =  size(output.paths{4},3);
%calculate the transition time
for i = 1:num_paths
time(i) =sum(output.paths{4}(:,1,i)<0)/M;
end
% plot the result
plot(1:num_paths,time)
title('Transition Times','FontSize',16)
axis([1,num_paths,0,1])

%make a plot of autocorrelation of path{4}
%compute the number of lag times we can compute with our Fourier method
nLags = (cond.samples-cond.burn)/2;
%initialize a figure to hold the autocorrelation
autocorr=figure(3);
set(gca, 'fontsize', 40)
%use autocorrelation_mult function to calculate the autocorrelation
correlates = autocorrelation_mult(output.paths{4});
%plot the result
plot(0:nLags-1,correlates(:,1:nLags))
title('Autocorrelation Function (all times)','FontSize',16)
axis([0,nLags-1,-1,1])
toc