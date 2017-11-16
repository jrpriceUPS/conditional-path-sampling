clear all; close all; clc
addpath('../Functions')

%%%%%%%%%%%
%Dimension%
%%%%%%%%%%%
n = 2;

%%%%%%%%%%%%%%%%
%SDE parameters%
%%%%%%%%%%%%%%%%

%degree of Brownian noise
SDE.noise = 1/2;


%%%%%%%%%%%%%
%Drift terms%
%%%%%%%%%%%%%

%depth of wells
k = 5;

%standard quartic SDE drift and derivative of drift
%f = -grad(U), where U is the potential energy landscape
%the input x is a row vector of length n
U     = @(x) -2*k*x(1)^2 + x(2)^2 + (x(1)^2)*(x(2)^2) + k*(x(1)^4) + x(2)^4;
f{1}  = @(x) [4*k*x(1) - 2*x(1)*(x(2)^2) - 4*k*(x(1)^3), -2*x(2) - 2*(x(1)^2)*x(2) - 4*(x(2)^3)];
df{1} = @(x) [4*k - 2*(x(2)^2) - 12*k*(x(1)^2), -4*x(1)*x(2);
              -4*x(1)*x(2), -2 - 2*(x(1)^2) - 12*(x(2)^2)]; % negative Hessian of U

drifts.f = f;
drifts.df = df;
drifts.U = U;
drifts.k = k;


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
cond.mean       = [-1 0; 1 0];        %top row is start point, bottom row is end point.
cond.cov_start  = [.0001 0; 0 .0001];   %(invertible) covariance matrices of start
cond.cov_end    = [.0001 0; 0 .0001];   %and end distributions
cond.std        = [.01 .01];
%0 correlation, standard deviation of .01 in each variable

%functions relating to the potentials of initial and final observations
%negative log of a multivariate normal distribution
%this part needs to be optimized
cond.start_neg_log    =  @(x) (x - cond.mean(1,:))*inv(cond.cov_start)*(x - cond.mean(1,:)).'/2;
cond.end_neg_log      =  @(x) (x - cond.mean(1,:))*inv(cond.cov_end)*(x - cond.mean(1,:)).'/2;
cond.start_d_neg_log  =  @(x) (inv(cond.cov_start) + inv(cond.cov_start)')*(x - cond.mean(1,:)).'/2;
cond.end_d_neg_log    =  @(x) (inv(cond.cov_end) + inv(cond.cov_end)')*(x - cond.mean(1,:)).'/2;

%draw an initial position from the initial distribution
cond.initial_pos = cond.mean(1,:) + cond.std(1)*randn(1,n);

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
HMC_params.L  = 300;    %set to 0 for it to be automatically calculated
HMC_params.rnd_L = 0; %amount of randomness in selecting HMC_params.L (percentage)

if ~HMC_params.L
    HMC_params.L  = round(5/HMC_params.dt);
end

%%%%%%%%%%%%%%%%%%%%%
%Plotting parameters%
%%%%%%%%%%%%%%%%%%%%%

%1 if plots should be generated, 0 otherwise
plots.show         =  0;   %1 if plots should be generated during simulation, 0 otherwise
plots.print_ratio  =  1;   %1 if acceptance rate should be printed after each step
plots.num_plotted  =  0;   %number of plots highlighted at end

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%Actual simulation%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%


output = conditional_path_mult(SDE,drifts,cond,domain,HMC_params,plots,n);
%output.accept_rate %prints accepte rate
