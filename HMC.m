function [out,accept] = HMC_mult(U,grad_U,dt_min,dt_max,L_min, L_max,current_q)
%
%Computes a new proposed vector using the hybrid Monte Carlo algorithm and
%accepts it or rejects it according to the Hamiltonian provided
%
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%U       =  the potential of the Hamiltonian
%dU      =  the gradient of the potential of the Hamiltonian
%dt_min  =  the smallest permissible timestep
%dt_max  =  the largest permissible timestep
%L_min   =  the smallest permissible number of timesteps to take
%L_max   =  the largest permissible number of timesteps to take
%q       =  the current vector from which to generate a new one
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%out     =  the new vector
%accept  =  a logical variable indicating whether the proposal was accepted
%           (1) or rejected (0)

%compute the resolution
M = size(current_q,1);

%compute the dimension
n = size(current_q,2);

%compute the timestep to use in this HMC step
dt = dt_min + (dt_max-dt_min)*rand;
L = round(L_min + (L_max - L_min)*rand);

%record the current vector
q = matvec(current_q); %"unraveling" the matrix (columnwise)into a column vector

%draw initial virtual momentum from the standard multivariate normal
%distribution
p = randn(size(q));
current_p = p;

%half step for momentum at beginning
q = vecmat(q,M,n); %converting q back to a matrix to be able to input it to grad_U
g1 = matvec(grad_U(q));
p = p - dt*g1/2;
q = matvec(q); %converting q back to a vector

%alternate full steps for position and momentum
for i=1:L

    %full step for positions
    q = q + dt * p;

    %make full step for the momentum except at the end of the trajectory
    if i~=L
       q = vecmat(q,M,n);
       g2 = matvec(grad_U(q));
       p = p - dt * g2; 
       q = matvec(q);
    end
    
end

%make a half step for momentum at the end
q = vecmat(q,M,n);
g3 = matvec(grad_U(q));
p = p - dt * g3/2;


%negate the momentum to make proposal symmetric
p = -p;

%evaluate potential and kinetic energies at start and end of trajectory
current_U = U(current_q);
current_K = sum(current_p.^2)/2;
proposed_U = U(q);
proposed_K = sum(p.^2)/2;

%accept or reject the state at the end based on this
if rand < exp(current_U-proposed_U+current_K-proposed_K)
    out = q; %accept
    accept = 1;
else
    out = current_q; %reject
    accept = 0;
end