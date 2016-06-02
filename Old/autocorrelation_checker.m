clear all;close all;

n  = 1;
dt = 0.01;
T = 1;
phi = 0.3;

t = 0:dt:T;



x=zeros(length(t),n);


for i=2:length(t)
    x(i,:) = phi*x(i-1,:) + sqrt(dt)*randn(1,n);
end

tic
corr1 = autocorrelation(x.');
toc

tic
corr2 = autocorrelation_v2(x.');
toc

[a,b]=size(corr1);

for i=1:n
figure(2*(i-1)+1)
plot(0:19,corr1(i,1:20))
hold on
plot(0:19,corr2(i,1:20),'k')
plot(0:19,phi.^(0:19),'r')

figure(2*i)
subplot(2,1,1)
plot(0:b-1,corr1(i,:)-phi.^(0:b-1));
subplot(2,1,2)
plot(0:b-1,corr2(i,:)-phi.^(0:b-1));
end