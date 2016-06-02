function out = gaussian_deriv(x,mu,sigma,w)
%
%A function to evaluate the derivative of a set of gaussians at positions x
%
%Each gaussian has the form:
%
%w * Exp[-(x - mu)^2 / (2 * sigma^2)]
%
%
%%%%%%%%
%Input:%
%%%%%%%%
%
%x      =  the location to evaluate (can be a vector)
%mu     =  the mean of each gaussian (a vector)
%sigma  =  the standard deviation of each gaussian (a vector)
%w      =  the weight of each gaussian
%
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%out  =  a vector of the derivative evaluated at the point(s) x

%compute the number of evaluation points
num_x     = length(x);

%compute the number of gaussians
num_gauss = length(mu);

%construct matrices to evaluate each gaussian at each position so they can
%be summed later
x_mat = repmat(x,1,num_gauss);
mu_mat = repmat(mu.',num_x,1);
sigma_mat = repmat(sigma.',num_x,1);
w_mat = repmat(w.',num_x,1);

%evaluate and sum the results
out = -sum(w_mat.*(x_mat-mu_mat)./sigma_mat.^2.*exp(-(x_mat-mu_mat).^2./(2*sigma_mat.^2)),2);