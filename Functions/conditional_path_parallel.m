function output = conditional_path_parallel(SDE,drifts,cond,domain,HMC_params,plots)
%
%A function to compute conditional paths of the SDE
%
%dX = f_i(X)dt + sigma(X)dB_t
%
%with initial distribituion g(X_0) and end distribution h(X_T), utilizing 
%path exchanges between adjacent drift levels f_i and f_{i+1} to 
%decorrelate samples
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
%drifts: a structure containing the modified and final drifts and related
%        parameters
%
%drifts.f   =  a cell array of function handles of the modified drifts
%drifts.df  =  a cell array of function handles of the derivatives of the
%              modified drifts
%
%
%cond: a structure containing initial and end distribution information and
%      modified drifts (the condition)
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
%output.paths  =  an array of conditional paths governed by the different 
%                 levels of drift
%                 Dimensions: (time) x (drift levels) x (independent samples)
%
%output.accept_rate_HMC   =  proportion of acceptances of HMC proposals at 
%                            each drift level
%
%output.accept_rate_swap  =  proportion of acceptances of swap proposals at
%                            each drift level

%load information about the stochastic differential equation
sigma  =  SDE.noise;

%load information about the modified drifts
f_list    =  drifts.f;
df_list   =  drifts.df;

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
dt  =  domain.dt;
T   =  domain.endtime;

%load information about the hybrid monte carlo step
dt_HMC  =  HMC_params.dt;
T_HMC   =  HMC_params.T;

%load information about plotting
show         =  plots.show;
subplot_dim  =  plots.subplot_dim;

%initialize the time domain
t=0:dt:T;

%initialize acceptance rates
accept_rate_HMC   =  zeros(length(f_list),1);
accept_rate_swap  =  zeros(length(f_list)-1,2);

%compute the number of HMC steps to take
L_HMC = T_HMC/dt_HMC;

%compute paths that should be kept (according to gap) and initialize
%counter
indices = burn+1:gap+1:samples;
k = 1;

%initialize the output containing the paths
paths = zeros(length(t),length(f_list),length(indices));

%compute initial paths if not provided
if isfield(SDE,'initial_path')
    
    initial_path  =  SDE.initial_path;
    
else
    
    %compute initial paths with Langevin dynamics
    initial_path       =  zeros(length(t),length(f_list));
    initial_path(1,:)  =  X0;
    
    for j=1:length(f_list)
        
        for i=1:T/dt
            
            initial_path(i+1,j) = initial_path(i,j) + f_list{j}(initial_path(i,j))*dt + sigma*sqrt(dt)*randn;
            
        end
        
    end
    
end

%initialize current path with the initial path
current_path = initial_path;

%create cell array of potentials and gradient potentials
U   =  cell(1,length(f_list));
dU  =  cell(1,length(f_list));

for i=1:length(f_list)
    
    U{i}   =  @(b) potential(b,start_obs_potential,end_obs_potential,f_list{i},sigma,dt);
    dU{i}  =  @(b) grad_potential(b,start_obs_d_potential,end_obs_d_potential,f_list{i},df_list{i},sigma,dt);
    
end

%if plotting, initialize figure with subplots
if show == 1
    
    figure(1)
    
    for i=1:length(f_list)
        subplot(subplot_dim(1),subplot_dim(2),i)
        title(sprintf('Level = %g',i));
    end
    
end



%take samples
for i=1:samples
  
    %use hybrid Monte Carlo to propose and accept or reject new paths for
    %each drift level
    for j=1:length(f_list)
        
        [newpath,accept]  =  HMC(U{j},dU{j},0.4*dt_HMC,1.2*dt_HMC,L_HMC,current_path(:,j));
        current_path(:,j) = newpath;
        accept_rate_HMC(j) = accept_rate_HMC(j) + accept/samples;
        
    end
    
    
    %randomly select a path to swap with its nearest neighbor
    swap = randi(length(f_list)-1,1);
    
    %record how many swaps were attempted at each level
    accept_rate_swap(swap,2) = accept_rate_swap(swap,2) + 1;
    
    %complete the swap with a Metropolis accept/reject step
    if rand < exp(U{swap}(current_path(:,swap))+U{swap+1}(current_path(:,swap+1))-U{swap}(current_path(:,swap+1))-U{swap+1}(current_path(:,swap))); 
        
        current_path(:,[swap,swap+1]) = current_path(:,[swap+1,swap]);
        
        %update the number of successes of this particular swap
        accept_rate_swap(swap,1) = accept_rate_swap(swap,1) + 1;
        
    end
    
    %record paths
    if i == indices(k)
        
        paths(:,:,k) = current_path;
        k = k+1;
        
        %avoid errors due to extra looping
        if k > length(indices)
            k = 1;
        end
        
        %if plotting, display paths for each drift level
        if show==1
            
            for j=1:length(f_list)
                
                subplot(subplot_dim(1),subplot_dim(2),j)
                hold all
                plot(t,current_path(:,j))
                
            end
            
            drawnow
            
        end
        
    end
    
end


%save results
output.paths             =  paths;
output.accept_rate_HMC   =  accept_rate_HMC;
output.accept_rate_swap  =  accept_rate_swap(:,1)./accept_rate_swap(:,2);