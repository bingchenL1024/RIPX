% Bingchen Liu May 21, 2025
% This code generate example forcing using SDE 

close all
clear all 

th = 0.1;
mu = 0;
sig = 0.3;
dt = 1e-2;
tmax = 20;
t = 0:dt:tmax;             % Time vector
x = zeros(1,length(t)); % Allocate output vector, set initial condition
rng(1);                 % Set random seed
for i = 1:length(t)-1
    x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end

[acf,lags] = autocorr(x);

figure()
subplot(121)
plot(t,x);
subplot(122)
plot(lags,acf)


%% Method in Jekins and Watts (using 0 mean of x)
close all
clear all 

tau = 1; % decorrelation time scale in autocrr of e^{-t/tau}
sig = 0.3;
dt = 1e-2;
tmax = 20;
t = 0:dt:tmax;             % Time vector
x = zeros(1,length(t)); % Allocate output vector, set initial condition
rng(1);                 % Set random seed
for i = 1:length(t)-1
    x(i+1) = x(i)-(dt/tau)*x(i)+(dt/tau)*randn;
end
x_mean =ones(length(t),1)*mean(x);

t_lag = 0:dt:tmax;
auto_pred = exp(-t_lag/tau);

[acf,lags] = autocorr(x);

%%
figure()
subplot(121)
plot(t,x);
hold on 
plot(t,x_mean,'r','LineWidth',5)
hold off 
subplot(122)
plot(lags,acf)
hold on 
plot(t_lag,auto_pred)
