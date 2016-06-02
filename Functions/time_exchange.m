function [x_new,y_new,accept] = time_exchange(x,y,sigma,U_x,U_y)
%
%Attempts a Metropolis path exchange between paths x and path y where path
%y is twice as long as path x:
%
%x = [u_2,u_4,...,u_{T-2},u_T]
%y = [u_1,u_2,...,u_{T-1},u_T]
%
%The length of y must be 2^l and the length of x must be 2^(l-1)

%timestep of the long path
dt_y = 1/(length(y)-1);

%proposed new short path is truncated version of old long path
x_hat = y(1:2:end);

%proposed new long path requires sampling for the in-between points 
r1 = randn(length(x)-1,1)*sigma*sqrt(dt_y) + 0.5*(x(1:end-1) + x(2:end));

%construct the proposed long path
y_hat = zeros(size(y));
y_hat(1:2:end)    =  x;
y_hat(2:2:end-1)  =  r1;

%metropolis algorithm requires sampling for the in-between points of the
%old long path for comparison
r2 = randn(length(x)-1,1)*sigma*sqrt(dt_y) + 0.5*(y(1:2:end-2) + y(3:2:end));

%construct the modified old long path
y_tilde             =  zeros(size(y));
y_tilde(1:2:end)    =  y(1:2:end);
y_tilde(2:2:end-1)  =  r2;

%log acceptance probability
prob = U_x(x)-U_x(x_hat)+U_y(y_tilde)-U_y(y_hat) ...
    + sum((r1 - 0.5*(x(2:end) + x(1:end-1))).^2/(2*sigma^2*dt_y)) ...
    - sum((r2 - 0.5*(y(3:2:end) + y(1:2:end-2))).^2/(2*sigma^2*dt_y));

if rand < exp(prob)
     
    %if success, exchange paths and mark success
    x_new = x_hat;
    y_new = y_hat;
    accept = 1;
    
else
    
    %if fail, leave paths as is and do not mark success
    x_new = x;
    y_new = y;
    accept = 0;
    
end