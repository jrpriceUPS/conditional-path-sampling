function out = gaussian_2nd_deriv(x,mu,sigma,w)

num_x     = length(x);
num_gauss = length(mu);

x_mat = repmat(x,1,num_gauss);

mu_mat = repmat(mu.',num_x,1);
sigma_mat = repmat(sigma.',num_x,1);
w_mat = repmat(w.',num_x,1);

out = sum(w_mat.*(mu_mat.^2-sigma_mat.^2-2*mu_mat.*x_mat+x_mat.^2)./(sigma_mat.^4).*exp(-(x_mat-mu_mat).^2./(2*sigma_mat.^2)),2);