function output = conditional_path_v2(SDE,drifts,cond,domain,HMC_params,plots)
%
%A function to compute conditional paths of the SDE
%
%dX = f(X)dt + sigma(X)dB_t
%
%with initial condition X_0 and end (noisy) observation X_T
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%SDE: a structure containing information about the stochastic differential
%equation
%
%SDE.noise         =  a function handle for the noise level sigma(X)
%SDE.initial       =  the initial condition X_0
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
%cond: a structure containing end distribution information and modified
%      drifts (the condition)
%
%cond.mean      =  the end observation mean X_T
%cond.var       =  the variance of the observation
%cond.samples   =  the number of samples to take at each level
%cond.burn      =  the number of samples at the start of the simulation to
%                  disregard
%cond.gap       =  how frequently to save results after the burn period
%                  (the number of unsaved paths between saved paths)
%
%
%domain: a structure detailing the domain of simulation
%
%domain.dt         =  the time step
%domain.endtime    =  the end time
%domain.periodic   =  a logical variable that is 1 if the spatial domain is
%                     periodic and 0 otherwise
%domain.endpoints  =  the left and right endpoints of a periodic domain
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
%plots.show        =  a logical variable that is 1 if plots should be displayed
%                     and 0 otherwise
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%output.paths  =  a matrix of paths governed by the original dynamics
%sampled from the conditional distribution


%load information about the stochastic differential equation
sigma  =  SDE.noise;
X0     =  SDE.initial;

%load information about the modified drifts
f_list   =  drifts.f;
df_list  =  drifts.df;

%load information about the conditional end point
XT         =  cond.mean;
obs_noise  =  cond.var;

%load information about sampling
samples  =  cond.samples;
gap      =  cond.gap;
burn     =  cond.burn;

%load information about the domain of simulation
dt        =  domain.dt;
T         =  domain.endtime;
periodic  =  domain.periodic;

%load information about the hybrid monte carlo step
dt_HMC  =  HMC_params.dt;
T_HMC   =  HMC_params.T;

%load information about plotting
show  =  plots.show;

%initialize the time domain
t=0:dt:T;

%initialize acceptance rate
accept_rate = zeros(length(f_list),1);

%compute the number of HMC steps to take
L_HMC = T_HMC/dt_HMC;

%initialize periodic parameters if needed
if periodic == 1
    L       =  domain.edges(1);
    R       =  domain.edges(2);
    period  =  R-L;
end

%compute indices of paths that should be independent (according to gap)
independent_indices = burn+1:gap+1:samples;
k = 1;

%initialize the output containing the paths
paths = zeros(length(t),length(independent_indices));
paths(1,:) = X0;

%compute an initial path if one is not provided
if isfield(SDE,'initial_path')
    
    current_path  =  SDE.initial_path;
    
else
    
    current_path     =  zeros(length(t),1);
    current_path(1)  =  X0;
    
    for i=1:T/dt
        current_path(i+1,end) = current_path(i,end) + f_list{1}(current_path(i,end))*dt + sigma*sqrt(dt)*randn;
        
        if periodic == 1
            if current_path(i+1,end)<L
                current_path(i+1,end) = current_path(i+1,end)+period;
            elseif current_path(i+1,end)>R
                current_path(i+1,end) = current_path(i+1,end)-period;
            end
        end
        
    end
    
end

%remove the (unchanging) initial condition
current_path = current_path(2:end);


for j=1:length(f_list);
    
    %update the modified drift and its derivative (in the last iteration,
    %this is the original drift)
    modified_drift = f_list{j};
    modified_d_drift = df_list{j};
    
    %update the potential energy and its derivative
    if periodic ~=1
        
        U = @(b) potential_v2(b,X0,XT,obs_noise,modified_drift,sigma,dt);
        dU = @(b) grad_potential_v2(b,X0,XT,obs_noise,modified_drift,sigma,modified_d_drift,dt);
        
    else
        
        U = @(b) periodic_potential_v2(b,X0,XT,obs_noise,modified_drift,sigma,dt,L,R);
        dU = @(b) periodic_grad_potential_v2(b,X0,XT,obs_noise,modified_drift,sigma,modified_d_drift,dt,L,R);
        
    end
    
    %take samples
    for i=1:samples
        
        [newpath,accept]  =  HMC(U,dU,0.4*dt_HMC,1.2*dt_HMC,L_HMC,current_path);
        current_path = newpath;
        accept_rate(j) = accept_rate(j) + accept/samples;
        
        if i == independent_indices(k)
            
            paths(2:end,k) = current_path;
            k = k+1;
            
        end
        
    end
    
    %plot if requested
    if show == 1
        
        figure(1)
        plot(t,paths)
        
        if j~=length(f_list)
            title(sprintf('Level %g',j))
        else
            title('Original Dynamics')
        end
        drawnow
        
    end
    
end

%save results
output.paths        =  paths;
output.accept_rate  =  accept_rate;