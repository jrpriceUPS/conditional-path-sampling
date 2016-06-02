function output = meta_processing(cond,SDE,data,meta)
%
%A function to post-process data from metadynamics simulation into a form
%useful for conditional path sampling, yielding cell arrays of a sequence
%of drifts (and their derivatives) from the end of the metadynamics to the
%original dynamics
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%cond: a structure containing end distribution information and modified
%      drifts (the condition)
%
%cond.n_levels  =  the number of levels of modified drifts to construct
%                  from the metadynamics data
%
%
%SDE: a structure containing information about the stochastic differential
%equation
%
%SDE.drift        =  the original drift term of the SDE
%SDE.drift_deriv  =  the derivative of the original drift term of the SDE
%
%
%data: a structure containing information about the metadynamics gaussian
%deposits
%
%data.locations  =  the locations of each gaussian deposit (in order)
%
%data. weights   =  the weights of each gaussian deposit (in order)
%
%
%meta: a structure containing metadynamical parameters
%
%meta.width  =  the width of each gaussian deposit
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%output.f   =  a cell array of n_levels of modified drifts and the original
%              drift
%output.df  =  a cell array of the derivatives of the above drifts


%load the number of levels to construct
n_levels = cond.n_levels;

%identify indices in order to evenly split the metadynamics biases into the
%n_levels
indexes = floor(length(data.locations)/n_levels):floor(length(data.locations)/n_levels):floor(length(data.locations));

%initialize output cells
f = cell(n_levels+1,1);
df = cell(n_levels+1,1);

%set the final drift and derivative to the original
f{n_levels+1} = SDE.drift;
df{n_levels+1} = SDE.drift_deriv;

%set the intermediate drifts and derivatives according to metadynamics
for i=1:n_levels
    this_index = n_levels+1-i;
    f{i}   =  @(y) f{n_levels+1}(y) - gaussian_deriv(y,data.locations(1:indexes(this_index)),meta.width*ones(indexes(this_index),1),data.weights(1:indexes(this_index)));
    df{i}  =  @(y) df{n_levels+1}(y) - gaussian_2nd_deriv(y,data.locations(1:indexes(this_index)),meta.width*ones(indexes(this_index),1),data.weights(1:indexes(this_index)));
end

%save results
output.f = f;
output.df = df;