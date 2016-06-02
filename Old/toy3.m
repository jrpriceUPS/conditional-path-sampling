close all;

%kB = 1.38064852*10^-23;
%dTemp = 10000;
kB = 1;
dTemp = 10000;


%timestep and endtime
dt = .1;
T = 150;

%width of gaussians
width = 0.3;

%initial gaussian size
w = 0.3;

%how frequently we deposit gaussians
howoften=10;

deposit_list = zeros(T/dt/howoften,1);
weight_list = zeros(T/dt/howoften,1);
index=1;



f_0 = @(x) 3*sin(x);
U_0 = @(x) 3*cos(x);

x=-pi;
f = f_0;
U = U_0;

U_mod = @(x) 0;

xgrid=linspace(-2*pi,2*pi,1000);

figure
for i=1:T/dt
    x = x + f(x)*dt + 1/2*sqrt(dt)*randn;
    
    if x<-2*pi
        x = 4*pi-x;
    elseif x>2*pi
        x = x-4*pi;
    end
    
    
    if mod(i,howoften) == 0
        deposit_list(index) = x;
        weight_list(index)  = w*exp(-U_mod(x)/(kB*dTemp));
        
        f_mod = @(y) -gaussian_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index))...
            -gaussian_deriv(y,deposit_list(1:index)+4*pi,repmat(width,1,index),weight_list(1:index))...
            -gaussian_deriv(y,deposit_list(1:index)-4*pi,repmat(width,1,index),weight_list(1:index));
        U_mod = @(y) gaussian(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index))...
            +gaussian(y,deposit_list(1:index)+4*pi,repmat(width,1,index),weight_list(1:index))...
            +gaussian(y,deposit_list(1:index)-4*pi,repmat(width,1,index),weight_list(1:index));
        
        
        f = @(y) f_0(y) + f_mod(y);
        U = @(y) U_0(y) + U_mod(y);
        
        
        index = index + 1;
    end
    
    plot(xgrid,U(xgrid),'r')
    hold on
    plot(x,U(x),'r.')
    plot(xgrid,U_0(xgrid))
    axis([-2*pi,2*pi,-4,10])
    hold off
    drawnow
end












samples=10;
n_modified=10;

X0 = -pi;
XT = pi;
obs_noise = 0.01;
T = 20;

dt = .1;
t=0:dt:T;

dt_HMC = 0.01;
T_HMC = .1/dt_HMC;

increments = zeros(length(t)-1,samples+1);
increments(:,end) = 1/2*randn(length(t)-1,1);


indexlist=0:ceil(length(deposit_list)/n_modified):length(deposit_list);

for j=1:length(indexlist);
    
    index=indexlist(length(indexlist)-j+1);
    
    increments(:,1) = increments(:,end);
    modified_drift = @(y) f_0(y)-gaussian_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index))...
        -gaussian_deriv(y,deposit_list(1:index)+4*pi,repmat(width,1,index),weight_list(1:index))...
        -gaussian_deriv(y,deposit_list(1:index)-4*pi,repmat(width,1,index),weight_list(1:index));
    modified_d_drift = @(y) 3*cos(y)...
            - gaussian_2nd_deriv(y,deposit_list(1:index),repmat(width,1,index),weight_list(1:index))...
            - gaussian_2nd_deriv(y,deposit_list(1:index)+4*pi,repmat(width,1,index),weight_list(1:index))...
            - gaussian_2nd_deriv(y,deposit_list(1:index)-4*pi,repmat(width,1,index),weight_list(1:index));
    
    U = @(b) periodic_potential(b,X0,XT,obs_noise,modified_drift,dt,-2*pi,2*pi);
    dU = @(b) periodic_grad_potential(b,X0,XT,obs_noise,modified_drift,modified_d_drift,dt,-2*pi,2*pi);
    
    for i=1:samples
        increments(:,i+1) = HMC(U,dU,dt_HMC,T_HMC,increments(:,i));
    end
    paths = zeros(length(t),samples+1);
    paths(1,:) = X0;


for k=1:samples+1
    for l=2:length(t)
        paths(l,k) = paths(l-1,k)+dt*modified_drift(paths(l-1,k))+increments(l-1,k);
        
        if paths(l,k) <-2*pi
            paths(l,k) = 4*pi+paths(l,k);
        elseif paths(l,k) > 2*pi
            paths(l,k) = paths(l,k)-4*pi;
        end
    end
end

figure(1)
plot(t,paths)
drawnow

end