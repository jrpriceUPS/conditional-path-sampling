close all

%kB = 1.38064852*10^-23;
%dTemp = 10000;
kB = 1;
dTemp = 1000;


%timestep and endtime
dt = 0.01;
T = 10;

%width of gaussians
width = 0.25;

%initial gaussian size
w = 0.25;

%how frequently we deposit gaussians
howoften=10;

deposit_list = zeros(T/dt/howoften,1);
weight_list = zeros(T/dt/howoften,1);
index=1;



f_0 = @(x) 5*gaussian_deriv(x,2/0.75,1/sqrt(2),1)+10*gaussian_deriv(x,-2/0.75,1/sqrt(2),1) - (4*(x+4).^3).*(x<-4.3193) - (4*(x-4).^3).*(x>4.25882);
U_0 = @(x) -5*gaussian(x,2/0.75,1/sqrt(2),1)-10*gaussian(x,-2/0.75,1/sqrt(2),1) + ((x+4).^4-1.690133).*(x<-4.3193) + ((x-4).^4-0.845067).*(x>4.25882);

x=-2/0.75;
f = f_0;
U = U_0;

U_mod = @(x) 0;

xgrid=linspace(-4,4,1000);

figure
for i=1:T/dt
    x = x + f(x)*dt + 1/2*sqrt(dt)*randn;
    
    
    if mod(i,howoften) == 0
        deposit_list(index) = x;
        weight_list(index)  = w*exp(-U_mod(x)/(kB*dTemp));
        
        f_mod = @(y) -gaussian_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index));
        U_mod = @(y) gaussian(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index));
        
        
        f = @(y) f_0(y) + f_mod(y);
        U = @(y) U_0(y) + U_mod(y);
        
        
        index = index + 1;
    end
    
    plot(xgrid,U(xgrid),'r')
    hold on
    plot(x,U(x),'r.')
    plot(xgrid,U_0(xgrid))
    axis([-4,4,-10,5])
    hold off
    drawnow
end







samples=10;
n_modified=10;

X0 = -2/0.75;
XT = 2/0.75;
obs_noise = 0.01;
T = 100;

dt = .1;
t=0:dt:T;

dt_HMC = 0.01;
T_HMC = .01/dt_HMC;

increments = zeros(length(t)-1,samples+1);
increments(:,end) = 1/2*randn(length(t)-1,1);


indexlist=0:ceil(length(deposit_list)/n_modified):length(deposit_list);

for j=1:length(indexlist);
    
    index=indexlist(length(indexlist)-j+1);
    
    increments(:,1) = increments(:,end);
    modified_drift = @(y) f_0(y)-gaussian_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index));
    modified_d_drift = @(y) -5*gaussian_2nd_deriv(x,2/0.75,1/sqrt(2),1)-10*gaussian_2nd_deriv(x,-2/0.75,1/sqrt(2),1)...
        - (12*(x+4).^2).*(x<-4.3193) - (12*(x-4).^2).*(x>4.25882)...
            - gaussian_2nd_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index));
    
    U = @(b) potential(b,X0,XT,obs_noise,modified_drift,dt);
    dU = @(b) grad_potential(b,X0,XT,obs_noise,modified_drift,modified_d_drift,dt);
    
    for i=1:samples
        increments(:,i+1) = HMC(U,dU,dt_HMC,T_HMC,increments(:,i));
    end
    paths = zeros(length(t),samples+1);
    paths(1,:) = X0;


for k=1:samples+1
    for l=2:length(t)
        paths(l,k) = paths(l-1,k)+dt*modified_drift(paths(l-1,k))+increments(l-1,k);
    end
end

figure(1)
plot(t,paths)
drawnow

end
