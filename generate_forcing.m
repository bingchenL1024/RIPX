% Bingchen Liu May 21, 2025
% This code generate example forcing using SDE 

close all
clear all 
%% Gille note 
clear
A=randn(1000,100);
for i=1:100
AcA(:,i)=xcov(A(:,i),A(:,i),'unbiased'); % autocovariance for
end
fAcA=fft(AcA(500:1500,:)); % Fourier transform of autocovariance
frequency=(0:500)/1000;
loglog(frequency,2*abs(mean(fAcA(1:501,:),2)),'LineWidth',3)
set(gca,'FontSize',16)
xlabel('Frequency (cycles per data point)','FontSize',16)
ylabel('Spectral energy','FontSize',16)


mean_AcA=mean(AcA,2);
fmean_AcA=fft(mean_AcA(500:1500));
hold on
loglog(frequency,2*abs(fmean_AcA(1:501,:))*1.1,'r','LineWidth',3)


a=A(:);
N=length(a);
aca=xcov(a,a,'unbiased');
faca=fft(aca(N-500:N+500));
loglog(frequency,2*abs(faca(1:501)),'LineWidth',3)
legend('average of FFTs of many autocovariances',...
'FFT of averaged autocovariance (scaled by 1.1)','without chunkify')

%% simple non-stochastic test
close all
clear all 


dt = 0.1;
tmax=  50;
t = 0:dt:tmax;
N = length(t);
N_half= (N-1)/2;
omega1  = 2*pi/10;
omega2  = 2*pi/5;
x = sin(t*omega1)+sin(t*omega2);



[covar,lag] = xcorr(x,'unbiased');
lag = lag*dt;
%x_fft= fft(x);
covar_fft= fft(covar(N_half:N_half+N-1)); % only take the middle section of the autocovar
spec_covar = 2*abs(covar_fft(1:(length(covar_fft)-1)/2+1))*dt; %multiple 2 for the other half of the fft

freq_fromlag = (0:(N-1)/2)./tmax; %

[freq,spec] = fft_data(x,dt);



figure()
subplot(221)
plot(t,x)
subplot(222)
plot(lag,covar)
subplot(223)
plot(freq,spec)
xlim([0,0.3])
subplot(224)
plot(freq_fromlag,spec_covar)
xlim([0,0.3])
%% Orstein-Uhlenbeck process 
clear
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

tau = 2; % decorrelation time scale in autocrr of e^{-t/tau}
sig = 0.3;
dt = 1e-3;
tmax = 200;
t = 0:dt:tmax;             % Time vector
N = length(t);
N_half= (N-1)/2;
x = zeros(1,length(t)); % Allocate output vector, set initial condition
%rng(1);                 % Set random seed
for i = 1:length(t)-1
    x(i+1) = x(i)-(dt/tau)*x(i)+(dt/tau)*randn(1,1);
end

xmean = mean(x);
x_mean_plot =ones(length(t),1)*xmean;
x= detrend(x); %demean the data


[covar,lag] = xcorr(x,'unbiased');
lag = lag*dt;
%x_fft= fft(x);
covar_fft= fft(covar(N_half:N_half+N-1)); % only take the middle section of the autocovar
spec_covar = 2*abs(covar_fft(1:(length(covar_fft)-1)/2+1))*dt; %multiple 2 for the other half of the fft
freq_fromlag = (0:(N-1)/2)./tmax; %

[freq,spec] = fft_data(x,dt);



%% compare with theory 
t_lag = 0:dt:tau*3;
auto_pred = exp(-t_lag/tau);
[acf,lags_num] = autocorr(x,'NumLags',2*(tau/dt));
lags = lags_num*dt;


%% plot 
linewd = 2;
figure()
subplot(121)
plot(t,x);
hold on 
yline(xmean,'LineWidth',linewd)
hold off 

subplot(122)
plot(lags,acf,'LineWidth',linewd)
hold on 
plot(t_lag,auto_pred,'LineWidth',linewd)
hold off 


figure()
subplot(221)
plot(t,x)
subplot(222)
plot(lag,covar)
subplot(223)
loglog(freq,spec)
xlim([0,0.3])
subplot(224)
loglog(freq_fromlag,spec_covar)
xlim([0,0.3])