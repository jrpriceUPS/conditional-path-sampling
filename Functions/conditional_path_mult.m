function output = conditional_path_mult(SDE,drifts,cond,domain,HMC_params,plots,n)

%load information about the stochastic differential equation
sigma  =  SDE.noise;

%load information about the drift
f    =  drifts.f;
df   =  drifts.df;
V    =  drifts.U; %free energy surface
k    =  drifts.k; %depth of wells; for plotting the free energy surface

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
dt  = domain.dt;
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

%initialize the time domain
t = 0:dt:T;
M = length(t); %resolution

%initialize acceptance rate
accept_rate = 0;

%compute paths to be kept (according to gap) and initialize counter
indices = burn+1:gap+1:samples;

%initialize the output containing the paths
paths = zeros(M,n,length(indices)); %paths(x,y,z) gives the yth coordinate of the xth point on the zth path

%compute an initial path if one is not provided
if isfield(SDE,'initial_path')
    
    initial_path  =  SDE.initial_path;
    
else

    %compute initial paths with Langevin dynamics
    initial_path       =  zeros(length(t),n); %ith column is the ith coordinate; each row is a point in R^n
    initial_path(1,:)  =  X0;

    for i = 1:(T/dt)

        initial_path(i+1,:) = initial_path(i,:) + f{1}(initial_path(i,:))*dt + sigma*sqrt(dt)*randn(1,n);

    end
end

%initialize current path with the initial path
current_path  =  initial_path;

%initialize potential and gradient of potential
U = @(b) potential_mult(b,start_obs_potential,end_obs_potential,f{1},sigma,dt);
dU = @(b) grad_potential_mult(b,start_obs_d_potential,end_obs_d_potential,f{1},df{1},sigma,dt);

%stores the acceptance probability at each step
percent_accepted = zeros(samples,1);

%take samples
for i=1:samples
    
    [newpath,accept]  =  HMC_mult(U,dU,(1-rnd_dt)*dt_HMC,(1+rnd_dt)*dt_HMC,(1-rnd_L)*L_HMC, (1+rnd_L)*L_HMC,current_path);
    current_path = newpath;
    accept_rate = accept_rate + accept/samples;
    percent_accepted(i) = accept_rate * samples / i;
    
    if print_ratio
        percent_accepted(i) %prints acceptance ratio at current step
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


% save percent_accepted.dat percent_accepted -ascii %saves the acceptance ratio at each step
% save paths.dat output.paths -ascii




% %plotting the free energy landscape
% [X,Y] = meshgrid(-2:.1:2);                                
% Z = @(a) -2*a*(X.^2) + (Y.^2) + (X.^2).*(Y.^2) + a*(X.^4) + Y.^4;
% surf(X,Y,Z(k))
% % view(2) %for a top down view
% 
% hold on
% 
% %plotting the initial path in 2D (projection onto the plane z = 0)
% % plot(initial_path(:,1)', initial_path(:,2)', 'r')
% % axis([-2,2,-2,2])
% 
% z = zeros(1,length(t)); %height of the path at each point
% for i = 1:length(t)
%     z(i) = V([initial_path(i,1), initial_path(i,2)]);
% end
% 
% %plotting the initial path in 3D
% plot3(initial_path(:,1)', initial_path(:,2)', z, 'r')
% axis([-2,2,-2,2,-1*k - 1,5,0,5]);