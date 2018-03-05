function output = conditional_path_time_mult(SDE,drifts,cond,domain,HMC_params,plots,n)
%path exchanges between paths of different resolutions (powers of two)

%load information about the stochastic differential equation
sigma  =  SDE.noise;

%load information about the drift
f    =  drifts.f;
df   =  drifts.df;
V    =  drifts.U; %free energy surface
%k    =  drifts.k; %depth of wells; for plotting the free energy surface

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
dt  = domain.dt;% now it's a vector
T   = domain.endtime;
X0  =  cond.initial_pos;

%load information about the hybrid monte carlo step
dt_HMC  =  HMC_params.dt;
L_HMC   =  HMC_params.L;
rnd_dt  =  HMC_params.rnd_dt;
rnd_L   =  HMC_params.rnd_L;

%load information about plotting
show  =  plots.show;
print_ratio = plots.print_ratio;
if show==1
    subplot_dim  =  plots.subplot_dim;
end;

%record the length(dt) for convience, which is the number of different resolutions 
num_sub = length(dt);

%initialize the time domain
t = cell(num_sub,1);

%compute the time domain of different resolution
for i=1:num_sub
    t{i}=0:dt(i):T;
end

%initialize acceptance rate
accept_rate_HMC   =  zeros(num_sub,1);
accept_rate_swap  =  zeros(num_sub-1,2);

%compute paths to be kept (according to gap) and initialize counter
indices = burn+1:gap+1:samples;
k = 1;

%initialize the output containing the paths
paths = cell(1,num_sub);

%create the paths of different lengths 
for i=1:num_sub    
    paths{i} = zeros(length(t{i}),n,length(indices));    
end

%compute an initial path if one is not provided
if isfield(SDE,'initial_path')    
    initial_path  =  SDE.initial_path;    
else

    %compute initial paths with Langevin dynamics
    initial_path     =  cell(1,num_sub);   

   for j=1:num_sub
        %compute initial paths with Langevin dynamics 
        initial_path{j} = zeros(length(t(j)),n);%ith column is the ith coordinate; each row is a point in R^n
        initial_path{j}(:,:) = X0(:,:,j);
        
        for i=1:T/dt(j)            
            initial_path{j}(i+1,:) = initial_path{j}(i,:) + 2*dt(j).*[1 0];
        end       
    end
end

%initialize current path with the initial path
current_path  =  initial_path;

U   =  cell(1,num_sub);
dU  =  cell(1,num_sub);

%initialize potential and gradient of potential
for i=1:num_sub
    U{i} = @(b) potential_mult(b,start_obs_potential,end_obs_potential,f{1},sigma,dt(i));
    dU{i} = @(b) grad_potential_mult(b,start_obs_d_potential,end_obs_d_potential,f{1},df{1},sigma,dt(i));
end

%stores the acceptance probability at each step
percent_accepted = zeros(samples,num_sub);

%if plotting, initialize figure with subplots
if show == 1    
    figure
    for i=1:num_sub
        subplot(subplot_dim(1),subplot_dim(2),i)
        title(sprintf('Resolution = 2^{%g}', round(log2(dt(i)))));
    end    
end

%make the landscape
[X,Y] = meshgrid(-2:.1:2);
Z = @(a) -2*a*(X.^2) + (Y.^2) + (X.^2).*(Y.^2) + a*(X.^4) + Y.^4;

for j=1:num_sub
    subplot(plots.subplot_dim(1),plots.subplot_dim(2),j)
    hold all
    %plotting the free energy landscape  
    surf(X,Y,Z(drifts.k))
    hold on
end

%take samples
for i=1:samples
    %use hybrid Monte Carlo to propose and accept or reject new paths for
    %each resolution
    for j=1:num_sub
        [newpath,accept]  =  HMC(U{j},dU{j},(1-rnd_dt)*dt_HMC(j),(1+rnd_dt)*dt_HMC(j),(1-rnd_L)*L_HMC(j), (1+rnd_L)*L_HMC(j),current_path{j});
        current_path{j} = newpath;
        accept_rate_HMC(j) = accept_rate_HMC(j) + accept/samples;
        percent_accepted(i,j) = accept_rate_HMC(j) * samples / i;
    
        if print_ratio
            percent_accepted(i,:) %prints acceptance ratio at current step
        end
    
        %randomly select a path to swap with its nearest neighbor
        swap = randi(num_sub-1,1);
    
        %record how many swaps were attempted at each level
        accept_rate_swap(swap,2) = accept_rate_swap(swap,2) + 1;
    
        %complete the swap with a Metropolis accept/reject step
        [x,y,accept] = time_exchange_mult(current_path{swap},current_path{swap+1},sigma,U{swap},U{swap+1},n);
        current_path{swap} = x;
        current_path{swap+1} = y;
    
        %update the number of successes of this particular swap
        accept_rate_swap(swap,1) = accept_rate_swap(swap,1) + accept;
    end
    
    %record paths
    if i == indices(k)     
        for n=1:num_sub        
            paths{n}(:,:,k) = current_path{n};
        end
        k = k+1;
        
        %avoid errors due to extra looping
        if k > length(indices)
            k = 1;
        end
        
        %if plotting, display paths for each drift level
        if show==1
            if k ~=1 
                for j=1:num_sub
                    subplot(subplot_dim(1),subplot_dim(2),j)
                    hold all
                    %plotting the initial path in 2D (projection onto the plane z = 0)
                    % plot(initial_path(:,1)', initial_path(:,2)', 'r')
                    % axis([-2,2,-2,2])
                    M = size(current_path{j},1); %resolution
                    height = zeros(1,M); %height of the path at each point
                    for ii = 1:M
                        height(ii) = U{j}(current_path{j}(ii,:));
                    end
                    
                    %plotting paths in 3D
                    plot3(current_path{j}(:,1), current_path{j}(:,2), height)
                    hold on
                end
            end
            drawnow
        end
    end
end
    
%save results
output.paths        =  paths;
output.accept_rate  =  accept_rate_HMC;