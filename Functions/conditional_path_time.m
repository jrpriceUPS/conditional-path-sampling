function output = conditional_path_time(SDE,drifts,cond,domain,HMC_params,plots)
%
%A function to compute conditional paths of the SDEs
%
%dX = f_i(X)dt + sigma(X)dB_t
%
%with initial distribituion g(X_0) and end distribution h(X_T), utilizing
%path exchanges between paths of different resolutions (powers of two)
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
%drifts.df2  =  the second derivative of the drift of the SDE
%
%
%cond: a structure containing start and end distribution information
%
%cond.implicit         =  logical variable indicating whether to use the
%                         implicit integration scheme (cond.implicit = 1)
%                         or the explicit scheme (cond.implicit = 0)
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
%domain.length(dt)  =  number of levels of resolution
%domain.dt          =  the time step of the highest resolution
%domain.endtime     =  the end time
%
%
%HMC: a structure containing hybrid Monte Carlo parameters
%
%HMC.dt  =  the time step in the HMC
%HMC.T   =  the end time of each HMC simulation
%
%
%plots: a structure detailing plotting parameters
%
%plots.show         =  a logical variable that is 1 if plots should be displayed
%                      and 0 otherwise
%plots.subplot_dim  =  the dimensions of subplots for plotting results
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%output.paths  =  a cell array of conditional paths governed by the different
%                 levels of drift (each cell is a different resolution)
%                 Dimensions: (time) x (independent samples)
%
%output.accept_rate_HMC   =  proportion of acceptances of HMC proposals at
%                            each resolution
%
%output.accept_rate_swap  =  proportion of acceptances of swap proposals at
%                            each resolution

%load information about the stochastic differential equation
sigma  =  SDE.noise;

%load information about the modified drifts
f    =  drifts.f;
df   =  drifts.df;
if cond.implicit==1
    df2  =  drifts.df2;
end

%load information about the conditional start and end points
start_obs_potential    =  cond.start_neg_log;
start_obs_d_potential  =  cond.start_d_neg_log;
X0                     =  cond.initial_pos;
end_obs_potential      =  cond.end_neg_log;
end_obs_d_potential    =  cond.end_d_neg_log;

%load information about sampling
samples  =  cond.samples;
gap      =  cond.gap;
burn     =  cond.burn;

%load information about the domain of simulation
T           =  domain.endtime;

%initialize timesteps (resolutions)
dt = domain.dt;

%load information about the hybrid monte carlo step
dt_HMC  =  HMC_params.dt;
L_HMC   =  HMC_params.L;

%load information about plotting
print_ratio  =  plots.print_ratio;
show         =  plots.show;
if show==1
    subplot_dim  =  plots.subplot_dim;
end;

%initialize the time domain
t = cell(length(dt),1);

for i=1:length(dt)
    
    t{i}=0:dt(i):T;
    
end

%initialize acceptance rates
accept_rate_HMC   =  zeros(length(dt),1);
accept_rate_swap  =  zeros(length(dt)-1,2);

%compute paths that should be kept (according to gap) and initialize
%counter
indices = burn+1:gap+1:samples;
k = 1;

%initialize the output containing the paths
paths = cell(1,length(dt));

for i=1:length(dt)
    
    paths{i} = zeros(length(t{i}),length(indices));
    
end

%compute initial paths if not provided
if isfield(SDE,'initial_path')
    
    initial_path  =  SDE.initial_path;
    
else
    
    initial_path     =  cell(length(dt),1);
    
    %compute initial paths with Langevin dynamics
    for j=1:length(dt)
        
        initial_path{j} = zeros(length(t{j}),1);
        initial_path{j}(1) = X0(j);
        
        for i=1:T/dt(j)
            
            initial_path{j}(i+1) = initial_path{j}(i) + f{1}(initial_path{j}(i))*dt(j) + sigma*sqrt(dt(j))*randn;
            
        end
        
    end
    
    
end

%initialize current path with the initial path
current_path = initial_path;

%create cell array of potentials and gradient potentials (using correct
%integration scheme)
U   =  cell(1,length(dt));
dU  =  cell(1,length(dt));

for i=1:length(dt)
    
    if cond.implicit==1
        
        U{i}   =  @(b) potential_implicit(b,start_obs_potential,end_obs_potential,f{1},df{1},sigma,dt(i));
        dU{i}  =  @(b) grad_potential_implicit(b,start_obs_d_potential,end_obs_d_potential,f{1},sigma,df{1},df2{1},dt(i));
        
    else
        
        U{i}   =  @(b) potential(b,start_obs_potential,end_obs_potential,f{1},sigma,dt(i));
        dU{i}  =  @(b) grad_potential(b,start_obs_d_potential,end_obs_d_potential,f{1},df{1},sigma,dt(i));
        
    end
    
end

%if plotting, initialize figure with subplots
if show == 1
    
    figure
    
    for i=1:length(dt)
        subplot(subplot_dim(1),subplot_dim(2),i)
        title(sprintf('Resolution = 2^{%g}', round(log2(dt(i)))));
    end
    
end

%stores the acceptance probability at each step
percent_accepted = zeros(samples, length(dt));

%take samples
for i=1:samples
    
    %use hybrid Monte Carlo to propose and accept or reject new paths for
    %each resolution
    for j=1:length(dt)
        
        [newpath,accept]  =  HMC(U{j},dU{j},.99*dt_HMC(j),1.01*dt_HMC(j),L_HMC(j),current_path{j});
        current_path{j} = newpath;
        accept_rate_HMC(j) = accept_rate_HMC(j) + accept/samples;
        percent_accepted(i,j) = accept_rate_HMC(j) * samples / i;
        
    end
    
    if print_ratio
        percent_accepted(i,:) %prints acceptance ratio at current step
    end
    
    
    %randomly select a path to swap with its nearest neighbor
    swap = randi(length(dt)-1,1);
    
    %record how many swaps were attempted at each level
    accept_rate_swap(swap,2) = accept_rate_swap(swap,2) + 1;
    
    %complete the swap with a Metropolis accept/reject step
    [x,y,accept] = time_exchange(current_path{swap},current_path{swap+1},sigma,U{swap},U{swap+1});
    current_path{swap} = x;
    current_path{swap+1} = y;
    
    %update the number of successes of this particular swap
    accept_rate_swap(swap,1) = accept_rate_swap(swap,1) + accept;
    
    
    %record paths
    if i == indices(k)
        
        for n=1:length(dt)
            paths{n}(:,k) = current_path{n};
        end
        k = k+1;
        
        %avoid errors due to extra looping
        if k > length(indices)
            k = 1;
        end
        
        %if plotting, display paths for each resolution
        if show==1
            
            if k ~=1
                
                for j=1:length(dt)
                    
                    subplot(subplot_dim(1),subplot_dim(2),j)
                    hold all
                    plot(t{j},paths{j}(:,k-1))
                    
                end
                
                drawnow
                
            end
            
        end
        
    end
    
end


%save results
output.paths             =  paths;
output.accept_rate_HMC   =  accept_rate_HMC;
output.accept_rate_swap  =  accept_rate_swap(:,1)./accept_rate_swap(:,2);

output.accept_rate_swap %prints the swap rate
%idea: highlight the swapped paths? Keep track of which paths they were?