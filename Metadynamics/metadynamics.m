function output = metadynamics(SDE,meta,domain,plots)
%
%A function to simulate metadynamics of the stochastic differential
%equation
%
%dX = f(X)dt + sigma(X)dB_t
%
%with initial condition X_0
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%SDE: a structure containing information about the stochastic differential
%equation
%
%SDE.drift      =  a function handle for the drift term f(X) in the SDE
%SDE.potential  =  a function handle for the potential U(X) of the drift
%SDE.noise      =  a function handle for the noise level sigma(X)
%SDE.initial    =  the initial condition
%
%
%meta: a structure containing metadynamical parameters
%
%meta.weight    =  the initial amplitude of gaussian deposits
%meta.width     =  the width of gaussian deposits
%meta.tempered  =  a logical variable that is 1 if we are using
%                  well-tempered metadyanmics and 0 otherwise
%meta.freq      =  the frequency of depositions (number of time steps)
%dTemp          =  the temperature of the simulation (a tuning parameter)
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
%plots: a structure detailing plotting parameters
%
%plots.show        =  a logical variable that is 1 if plots should be displayed
%                     and 0 otherwise
%plots.axes        =  [xmin,xmax,ymin,ymax]
%plots.resolution  =  number of grid points to use when plotting potential
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%output.path       =  a vector of path taken by the particle during the 
%                     simulation
%output.locations  =  a vector of the locations at which gaussians were
%                     deposited
%output.weights    =  a vector of the weights assigned to each gaussian if
%                     well-tempered metadynamics were used


%physical parameter: Boltzmann's constant
kB = 1.38064852*10^-23;


%load information about the stochastic differential equation
f_0    =  SDE.drift;
U_0    =  SDE.potential;
sigma  =  SDE.noise;
X0     =  SDE.initial;

%load information about the metadynamics
w         =  meta.weight;
width     =  meta.width;
tempered  =  meta.tempered;
freq      =  meta.freq;
dTemp     =  meta.temp;

%load information about the domain of simulation
dt        =  domain.dt;
T         =  domain.endtime;

%load information about plotting
show        =  plots.show;
plot_axes   =  plots.axes;
resolution  =  plots.resolution;



%initialize modifed potential and weight list if using well-tempered
%metadynamics
if tempered==1
    U_mod = @(x) 0;
    weight_list  =  zeros(T/dt/freq,1);
end

%initialize plotting data if plotting
if show==1
    figure
    xgrid = linspace(plot_axes(1),plot_axes(2),resolution);
    xgrid = xgrid.';
end

%initialize particle
x     =  zeros(T/dt+1,1);
x(1)  =  X0;

%initialize other output variables
deposit_list  =  zeros(T/dt/freq,1);
meta_index    =  1;

%initialize drift and potential
f  =  f_0;
U  =  U_0;


for i=1:T/dt
    x(i+1) = x(i) + f(x(i))*dt + sigma*sqrt(dt)*randn;
    
    if mod(i,freq) == 0
        deposit_list(meta_index) = x(i+1);
        
        if tempered==1
            weight_list(meta_index)  =  w*exp(-U_mod(x(i+1))/(kB*dTemp));
        else
            weight_list(meta_index)  =  w;
        end

            
        f_mod = @(y) -gaussian_deriv(y,deposit_list(1:meta_index),width*ones(meta_index,1),weight_list(1:meta_index));
        U_mod = @(y) gaussian(y,deposit_list(1:meta_index),width*ones(meta_index,1),weight_list(1:meta_index));
        
        
        f = @(y) f_0(y) + f_mod(y);
        U = @(y) U_0(y) + U_mod(y);
        
        
        meta_index = meta_index + 1;
    end
    
    if show==1
        plot(xgrid,U(xgrid),'r')
        hold on
        plot(x(i+1),U(x(i+1)),'r.')
        plot(xgrid,U_0(xgrid))
        axis(plot_axes)
        title(sprintf('T=%g',i*dt))
        drawnow
        hold off
    end
end

output.path       =  x;
output.locations  =  deposit_list;

if tempered==1
    output.weights  =  weight_list;
end