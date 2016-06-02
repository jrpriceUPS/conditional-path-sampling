close all;
samples=10;
n_modified=10;

X0 = -1;
XT = 1;
obs_noise = 0.001;
T = 1;

dt = 0.001;
t=0:dt:T;

dt_HMC = 0.001;
T_HMC = 0.01/dt_HMC;




xgrid=linspace(-2,2,1000);

d_drift = 1/n_modified;

drift_mod_list = d_drift:d_drift:1;


increments = zeros(length(t)-1,samples+1);
increments(:,end) = 1/2*randn(length(t)-1,1);

for j=1:length(drift_mod_list);
    
    increments(:,1) = increments(:,end);
    modified_drift = @(x) drift_mod_list(j)*drift(x);
    modified_d_drift = @(x) drift_mod_list(j)*d_drift_dt(x);
    
    for i=1:samples
        U = @(b) potential(b,X0,XT,obs_noise,modified_drift,dt);
        dU = @(b) grad_potential(b,X0,XT,obs_noise,modified_drift,modified_d_drift,dt);
        
        increments(:,i+1) = HMC(U,dU,dt_HMC,T_HMC,increments(:,i));
    end
    paths = zeros(length(t),samples+1);
    paths(1,:) = -1;

%compute the path
for k=1:samples+1
    for l=2:length(t)
        paths(l,k) = paths(l-1,k)+dt*modified_drift(paths(l-1,k))+increments(l-1,k);
    end
end

figure(1)
plot(t,paths)
drawnow

end


