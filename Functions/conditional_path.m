function output = conditional_path(SDE,drifts,cond,domain,HMC_params,plots)
%
%A function to compute conditional paths of the SDE
%
%dX = f(X)dt + sigma(X)dB_t
%
%with initial distribituion g(X_0) and end distribution h(X_T)
%
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%SDE: a structure containing information about the stochastic differential
%equation
%
%SDE.noise         =  a function handle for the noise level sigma(X)
%SDE.initial_path  =  (optional) an initial guess at the path
%
%
%drifts: a structure containing the drifts of the SDE and related parameters
%
%drifts.f    =  the drift of the SDE
%drifts.df   =  the derivative of the drift of the SDE
%
%
%cond: a structure containing start and end distribution information
%
%cond.start_d_neg_log  =  negative logarithm of the starting distribution 
%                         (function handle)
%cond.start_d_neg_log  =  the derivative of the above (function handle)
%cond.initial_pos      =  a vector of initial conditions drawn a priori 
%                         from the initial distribution (one for each drift 
%                         level)
%cond.end_neg_log      =  negative logarithm of the ending distribution 
%                         (function handle)
%cond.end_d_neg_log    =  the derivative of the above (function handle)
%cond.samples          =  the number of samples to take at each level
%cond.gap              =  how frequently to save results after the burn
%                         period (the number of unsaved paths between saved 
%                         paths, to save all set to zero)
%cond.burn             =  the number of samples at the start of the
%                         simulation to disregard
%
%
%domain: a structure detailing the domain of simulation
%
%domain.dt         =  the time step
%domain.endtime    =  the end time
%
%
%HMC_params: a structure containing hybrid Monte Carlo parameters
%
%HMC_params.dt  =  the time step in the HMC
%HMC_params.L   =  the number of time steps in the HMC
%
%
%plots: a structure detailing plotting parameters
%
%plots.show        =  a logical variable that is 1 if plots should be displayed
%                     and 0 otherwise
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%output.paths  =  a matrix of paths governed by the original dynamics
%                 sampled from the conditional distribution
%                 Dimensions: (time) x (independent samples)

%load information about the stochastic differential equation
sigma  =  SDE.noise;

%load information about the drift
f    =  drifts.f;
df   =  drifts.df;

%load information about the conditional start and end points
start_obs_potential    =  cond.start_neg_log;
start_obs_d_potential  =  cond.start_d_neg_log;
end_obs_potential      =  cond.end_neg_log;
end_obs_d_potential    =  cond.end_d_neg_log;

%load information about sampling
samples  =  cond.samples;
gap      =  cond.gap;
burn     =  cond.burn;

%load information about the domain of simulation
dt  =  domain.dt;
T   =  domain.endtime;
X0  =  cond.initial_pos;

%load information about the hybrid monte carlo step
dt_HMC  =  HMC_params.dt;
L_HMC   =  HMC_params.L;

%load information about plotting
show  =  plots.show;
print_ratio = plots.print_ratio;

%initialize the time domain
t=0:dt:T;

%initialize acceptance rate
accept_rate = 0;

%compute paths that should be kept (according to gap) and initialize
%counter
indices = burn+1:gap+1:samples;
k = 1;

%initialize the output containing the paths
paths = zeros(length(t),length(indices));

%compute an initial path if one is not provided
if isfield(SDE,'initial_path')
    
    initial_path  =  SDE.initial_path;
    
else
    
    %compute initial paths with Langevin dynamics
    initial_path     =  zeros(length(t),1);
    initial_path(1)  =  X0;
    
    for i=1:T/dt
        
        initial_path(i+1,end) = initial_path(i,end) + f{1}(initial_path(i,end))*dt + sigma*sqrt(dt)*randn;
        
    end
    
end

%initialize current path with the initial path
current_path  =  initial_path;

%initialize potential and gradient of potential
U = @(b) potential(b,start_obs_potential,end_obs_potential,f{1},sigma,dt);
dU = @(b) grad_potential(b,start_obs_d_potential,end_obs_d_potential,f{1},df{1},sigma,dt);

%stores the acceptance probability at each step
percent_accepted = zeros(samples,1);

%take samples
for i=1:samples
    
    [newpath,accept]  =  HMC(U,dU,.99*dt_HMC,1.01*dt_HMC,L_HMC,current_path);
    current_path = newpath;
    accept_rate = accept_rate + accept/samples;
    percent_accepted(i) = accept_rate * samples / i;
    
    if print_ratio
        percent_accepted(i) %prints acceptance ratio at current step
    end
    
    %record paths
    if i == indices(k)
        
        paths(:,k) = current_path;
        k = k+1;
        
        %avoid errors due to extra looping
        if k > length(indices)
            k = 1;
        end
        
        %if plotting, display paths for each drift level
        if show==1
            
            figure(1)
            hold all
            plot(t,current_path)
            drawnow
            
        end
        
    end
    
end



%save results
output.paths        =  paths;
output.accept_rate  =  accept_rate;

save percent_accepted.dat percent_accepted -ascii %saves the acceptance ratio at each step
% save paths.dat output.paths -ascii