Conditional Path Sampling README

There are four folders and a PDF in this directory containing everything needed to run conditional path sampling algorithms in Matlab. The PDF is an informal write-up of the current state of their project, written last May. I recommend reading through that to get a feel for our main goals and current practices.

In summary, we are looking to create new sampling methods to sample transition paths between two different locally stable states. In a real life example, this might be an unfolded protein and a folded protein. Both are “stable” in that they will not change rapidly from their current state. But there are transitions between these two states as the protein folds itself. These transition paths are of great interest to biochemical modelers.

In this project, we are working with a much simpler case, though it should be eventually generalizable given a larger computer. We work on a one dimensional, double-well potential here. We’d like to sample paths that begin in one well and pass into the other one. You can think about it as a particle sitting in this potential well that gets jostled every time step as if we are wiggling the potential.

We want to find physically possible paths that the particle could take from the bottom of one well into the other one as a result of those jostles. This will definitely happen if we wait long enough, but who wants to wait (especially since it could be years…)? We want to find ways to accelerate this path finding.

In the four folders are found the code base for this project.



FOLDER 1: Demos

This folder contains demos of different path sampling methods along with their analyses. Everything is well commented, so hopefully you can start to piece together how they work. These demos are things you want to run to actually generate results and explore how things change as you modify parameters.

basic_sampling_demo: this should be where you start. It is the simplest sampling algorithm you can come up with. The output will be a plot of sampled paths, a plot of the transition times, and a plot of the autocorrelation function. Our goal is to make the autocorrelation function decay quickly, and for the transition times to rapidly explore different transition times.

time_exchange_demo: an improvement of the algorithm proposed by Jonathan Weare. In this case, we simultaneously sample paths at different time discretizations, and use the metropolis accept-reject algorithm (from statistics) to exchange paths between discretizations and accelerate mixing.

path_exchange_demo: a different improvement in which we simultaneously sample paths with different (less steep) wells. The metropolis accept-reject algorithm is used to exchange paths.

metadynamics_exchange_demo: similar to the path_exchange_demo, except the different wells come from metadynamics. This would be the one that could most easily be generalized to a real-world application.



FOLDER 2: Metadynamics

This folder contains the functions needed to run the metadynamics algorithm (see the literature referenced in the summary PDF).

metadynamics: the main metadynamics function used for exploring a potential energy landscape by “filling the wells with sand as you go.”

meta_processing: a post-processing script to convert the data from metadynamics into a continuum of potential energy landscapes from the real landscape up to the “filled in” one.

gaussian: a function to evaluate a gaussian function with mean mu, standard deviation sigma, and weight w at points x.

gaussian_deriv: the x-derivative of the above gaussian.

gaussian_2nd_deriv: the second x-derivative of the above gaussian.



FOLDER 3: Functions

This folder contains all the functions associated with actually running the path-sampling algorithm.

conditional_path: main function for sampling conditional paths. Uses hybrid monte carlo to generate new paths from the current path, and saves the results as requested.

conditional_path_time: samples conditional paths using the time-exchange technique described by Jonathan Weare.

conditional_path_parallel: samples conditional path using path-exchange techniques between paths in different potential energy landscapes.

potential: computes the potential of a particular path (the negative log of the probability).

grad_potential: computes the gradient of the potential of a particular path.

potential_implicit: computes the potential of a particular path using implicit algorithm (have gotten strange results with this).

grad_potential_implicit: computes the gradient of the potential of a particular path using implicit algorithm (have gotten strange results with this).

time_exchange: metropolis exchange algorithm for two paths, one of which is twice as fine as the other. Used in the time-exchange algorithm.

HMC: hybrid monte carlo algorithm to generate a new path from the current path. Constructs a Hamiltonian system that should have the desired distribution as the steady state, and evolves the current state a set number of time steps.

autocorrelation: computes the autocorrelation function of a set of paths. Uses FFT to do so in an efficient way.



FOLDER 4: Old

Contains old files. I’ve kept it for posterity, but you shouldn’t need to use anything in here.