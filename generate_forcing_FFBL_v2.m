% Bingchen Liu July 2,2025
% Version: include slope bottom 
% This code generate the vorticity forcing using Falk's method with
% cross-shore: uniform forcing  
% alongshore: Fourier components in alongshore
% Time: AR1 in time in 3rd dimension

% Forcing decomposition G = Psi0*Psim1*Psim2 
% Psi0: background forcing with given Swb
% G0: magnitude from G0
% W: width/envolope function


addpath(genpath('/home/n2kumar/MAT_Mytoolbox/jlab/'));


clearvars -except test
close all

load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/data_example_cFbr_stats')

%% setup the simulation 
ind= 14; % use the example in snapshot simulation in paper fig 
SS = S(ind);
xb_run14 = SS.xb;
slp =0.02; %beach slope 


Ly=1500;
dy=1;
Lsz=xb_run14;
dx=1;
t_tot= 120;
dt=0.1;
x_max = -10;%shore line location so that the domain does go to 0 to prevent h=0
x=-Lsz:dx:x_max;
x_max_full = -259; %full x domain to mimic model data
x_full =x_max_full:dx:0 ;
nx=length(x);
nx_full = length(x_full);
y=0:dy:Ly;
ny=length(y);
t=0:dt:t_tot;
nt = length(t);
h_br = slp*Lsz;
h=x.*(-slp);
g= 9.81; %gravitational constant 
Tp = 8;
Hsb = SS.Hs_b;


%% parameters from given model simulation
x_run14 = SS.X2;
Dw_run14 = SS.dECG;
Dw = interp1(x_run14,Dw_run14,x);
Irb =slp/(SS.stpb)^0.5;
dirspr_b= SS.sigma_b;
ky0_nond = 9.5*10^(-4)*dirspr_b+1.9*10^(-3); % see paper equation
ky0 = ky0_nond./(Irb*h);
ky0_max =0.146; %cutoff for max ky0 allowed. Note ky0 blows up @h=0
ky0(ky0>ky0_max)  = ky0_max;
G0_nond= 0.041*dirspr_b+0.059; % see paper equation
G0 = G0_nond*(Dw./(h_br^2*(g.*h).^(1/2))).*(1/slp);
G0 = repmat(G0,ny,1,nt);

%%%%%%%%%%%%%%%%%%%%%%%%% original version (override for test)
% ky0_nond= 1.4*10^(-2);
% Irb=0.2;
% ky0 = ky0_nond./(Irb*h);


%%%%%%%%%%%%%%%%%%%%%%%%% QC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(121)
plot(x,ky0)
xlabel('x(m)')
ylabel('ky0')
subplot(122)
plot(x,G0(1,:,1))
xlabel('x(m)')
ylabel('G0')
%% generate ky
[X,Y]=meshgrid(x,y);

if mod(ny,2)==0
    ky=[-ny/2+1:1:ny/2 ]';    
else
    ky=[-(ny-1)/2:1:(ny-1)/2 ]';
end
ky=ky/Ly;
%ky=ky(abs(ky)<=0.2); %add ky limit to 0.2cmp
dky=ky(2)-ky(1);


%% get the magnitude function G0 first
mxx0=.075;
x00=100;
bb=x00/2^(1/3);
aa=mxx0*2^(1/3)*1.5/bb;
G0_matt=aa*(-X+x_max)./(1+((-X+x_max)/bb).^3); %G0 magnitude function
G0_matt = repmat(G0_matt,1,1,nt);
% plot(x,G0_matt(1,:,1))

%% ========> generate decorrelation time scale tau and c 

G0_ratio = G0(1,:,1)/max(G0(1,:,1));
tau = sqrt(Hsb/g) *(4.6*10^(-3)*exp(5.61*G0_ratio)+1.9); %constant tau parameterization 
%tau = 0.75*ones(length(x),1);tau=tau'; %constant tau parameterization 
%tau = sqrt(h/g);
c_phase = 1.1*sqrt(g.*h);
ctau = tau.*c_phase; %c*tau scaling for dx


%%%%%%%%%%%%%%%%%%%%%%%%% original version (override for test)
% G0_nond = G0_matt(1,:,1)/max(G0_matt(1,:,1));
% tau = sqrt(h/g).*(0.027*exp(4.37*G0_nond)+1.28);
% tau = 0.75*ones(length(x),1);tau=tau';
% %tau = sqrt(h/g);
% c_phase = sqrt(g.*h);
% ctau = tau.*c_phase; %c*tau scaling for dx

%% get Weibull spectra (S_wb)
%ky0=.02;%peF_comp alongshore wavenumber

phi=2*pi*rand(size(ky)); %random phase for the Fourier component(Note: needs to be the same for all x) 

for xloc = 1:nx
    S_wb(:,xloc)=make_WB(ky,ky0(xloc));
    var0=sum(S_wb(:,xloc))*dky; % the variance I think I should get
    S_wb(:,xloc) = S_wb(:,xloc) / var0; % the variance should now be 1 % becasue int (S_wb dky) = 1 here
    F_comp(:,xloc)=(2*ny*S_wb(:,xloc)).^.5.*exp(1i*phi); % Fourier component with random phase; Ny is for Matlab fft scaling
    Psi(:,xloc) = real(ifft(ifftshift(F_comp(:,xloc))));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%% QC %%%%%%%%%%%%%%%%%%%%%%%%
for xind = 1:nx
    Psi_var(xind) = var(Psi(:,xind));
end 

figure()
subplot(121)
pcolorcen(X,Y,Psi);
xlabel('x(m)')
ylabel('y(m)')
title('IFFT of Fourier compoments')
subplot(122)
plot(x,Psi_var)
xlabel('x(m)')
ylabel('var(Psi)')




%%%%%%%%%%%%%% check Psi using S_wb --> should be the same 
xloc= 150;
[f1,S1] = fft_data(Psi(:,1),dy);
[f2,S2] = fft_data(Psi(:,xloc),dy);
figure()
subplot(121)
loglog(f1,S1,ky,2*S_wb(:,1))
legend('loc1','Theory')
xlabel('k_y(m^{-1})')
ylabel('Spectra energy')
subplot(122)
loglog(f2,S2,ky,2*S_wb(:,xloc))
legend('loc1','Theory')
xlabel('k_y(m^{-1})')
ylabel('Spectra energy')


%% follow the wave crest and update all var 
nwave= round(t_tot/Tp);%total number waves want to include in the simulation. Generally >t_tot/Tp
active_wave = zeros(nwave,nt);
x_wave = -Lsz*ones(nwave,nt);
b_wave = zeros(nwave,nt); % b term that represents stochastic processes



for ind_wave = 1:nwave %initilize all waves only when it's time (set beginning properties)
    active_wave(ind_wave,round((ind_wave-1)*Tp/dt)+1:end)=1;
    b_wave(ind_wave,round((ind_wave-1)*Tp/dt)+1)=randn(1,1); %initilize b term for only the first step
end 


x_wave(:,1)=-Lsz; %initial locaiton of the wave crest 
for ind_wave = 1:nwave
    for tind = 1:nt-1
        if active_wave(ind_wave,tind) == 1 %only propagate wave when it's its turn
            x_wave(ind_wave,tind+1) = x_wave(ind_wave,tind)+ 1.1*dt*sqrt(g*(-slp*x_wave(ind_wave,tind))); %1.1 is the scaling factor
            if x_wave(ind_wave,tind+1) >x_max %if wave is in the domain 
            active_wave(ind_wave,tind+1:end) = 0; %deactive the wave if it's outside of the domain
            end 
        end 
    end 
end 

x_wave= round(x_wave);

% figure()
% imagesc(active_wave)
% xlabel('time (s)')
% ylabel('index for waves')
% title('Active map')
% 
% 
% figure()
% imagesc(x_wave)
% xlabel('time (s)')
% ylabel('index for waves')
% title('Wave location')
% cb=colorbar;
% cb.Label.String = 'Wave location (m)';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make b term %%%%%%%%%%%%%%%%%%%%%%%%%%
for ind_wave = 1:nwave
    for tind = 1:nt-1
        if active_wave(ind_wave,tind)==1 %only start propagate b forward if current b is known
            xloc = x_wave(ind_wave,tind);
            tau_loc = interp1(x,tau,xloc);
            b_wave(ind_wave,tind+1) = (1-dt/tau_loc)*b_wave(ind_wave,tind)+sqrt(2*dt/tau_loc)*randn(1,1);
        end 

    end 
end 




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% W (width function)%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = 3.16; %note need to divide by ctau in the loop for normalization
A2 = 3.57;
A3 = 1.82;
A4 = 1.1;
W_tot = zeros(ny,nx,nt);
W_wave_only_tot = zeros(ny,nx,nt);

for tind = 1:nt %@each time step 
    for ind_wave = 1:nwave
        if active_wave(ind_wave,tind) ==1
           x_wave_snap = x_wave(ind_wave,tind); % local wave location for given wave
           x_wave_snap=round(x_wave_snap);
           ctau_loc = interp1(x,ctau,x_wave_snap); %local ctau scaling 
           c_loc = interp1(x,c_phase,x_wave_snap); %local ctau scaling 
           lambda_crest = Tp*c_loc;%separation between wave crest 1.1 is the scaling factor in eq in paper
           %test_amp = (A1/(cos(A4)*sqrt(ctau_loc)));
           %W_wave = (A1/(cos(A4)*sqrt(ctau_loc))).*exp(-A2.*abs(X-x_wave_snap)/ctau_loc).*cos(A3.*(abs(X-x_wave_snap)/ctau_loc)+A4); %testing 
           W_wave = b_wave(ind_wave,tind)*(A1*sqrt(lambda_crest)/(cos(A4)*sqrt(ctau_loc))).*exp(-A2.*abs(X-x_wave_snap)/ctau_loc).*cos(A3.*abs(X-x_wave_snap)/ctau_loc+A4);%%correct one 
           W_tot(:,:,tind) = W_tot(:,:,tind)+W_wave;
           W_wave_only = (A1/(cos(A4)*sqrt(ctau_loc))).*exp(-A2.*abs(X-x_wave_snap)/ctau_loc).*cos(A3.*abs(X-x_wave_snap)/ctau_loc+A4);
           W_wave_only_tot(:,:,tind) = W_wave_only_tot(:,:,tind)+W_wave_only;
        end 
    end
    cFbr_temp(:,:,tind) = squeeze(W_tot(:,:,tind)).*Psi;
end 
W_wave_only_1s = W_wave_only_tot(:,:,1:1/dt:end); %subsample ;
W_1s = W_tot(:,:,1:1/dt:end); %subsample 
dim_w1s = size(W_1s);


%%%%%%%%%%%%%%%%%%%%%%%%%% QC %%%%%%%%%%%%%%%%%%%
% figure()
% for tind = 1:dim_w1s(3)
%     subplot(121)
%     plot(x,W_wave_only_1s(1,:,tind))
%     ylim([-3,3])
%     ylabel('W(x-ct)')
%     xlabel('x(m)')
%     subplot(122)
%     plot(x,W_1s(1,:,tind))
%     ylabel('W*b')
%     xlabel('x(m)')
%     ylim([-5,5])
%     drawnow
%     pause(0.5)
% end 


%% get final results (Finally!!!!)
cFbr_tot=G0.*cFbr_temp; %final field.
cFbr_sz = cFbr_tot(:,:,1:1/dt:end); %subsample it to 1s resolution. include only sz 
dim_cFbr_sz = size(cFbr_sz);
cFbr= zeros(ny,nx_full,dim_cFbr_sz(3));
cFbr(:,find(x_full ==x(1)):find(x_full==x(end)),:) = cFbr_sz; %put in the cFbr in the surfzone into the total matrix

[X_full,Y_full]=meshgrid(x_full,y);



%% get data stats 

dim_para= size(cFbr_sz);
for xind = 1:dim_para(2)
    for tind = 1:dim_para(3)
        cFbr_para_var(xind,tind) = sqrt(var(cFbr_sz(:,xind,tind)));
    end 
    cFbr_para_stats(xind) = mean(cFbr_para_var(xind,:));
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%% QC(G0 var)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(211)
plot(-258:0,cFbr_data_stats)
ylabel('mean_{t}(var_{y}(cFbr_{data}))')
xlim([-172,0])
subplot(212)
plot(x,cFbr_para_stats)
xlabel('x(m)')
ylabel('mean_{t}(var_{y}(cFbr_{para}))')
%% test plot for paperfig 


curlF_tstack = squeeze(cFbr(13,:,30:end));
curlF_snap = squeeze(cFbr(:,:,65));

dim_Hov = size(curlF_tstack);
x_Hov = -(dim_Hov(1)-1):0;
t_Hov = 0:dim_Hov(2)-1; %in seconds
[t_grid_Hov,x_grid_Hov] = meshgrid(t_Hov,x_Hov);
clim = 0.3;
figure()
subplot(211)
pcolorcen(X_full,Y_full,curlF_snap);
ylim([0,300])
cmocean('balance');
caxis([-clim,clim])
subplot(212)
pcolorcen(t_grid_Hov,x_grid_Hov,curlF_tstack);    
cmocean('balance');
caxis([-clim,clim])



%% save data 
head = 'data generated from generate_forcing_FFBL_V2, include realistic Dw from funwaveC';
save('/data1/bliu/data/parameterization_example_realistic','cFbr','curlF_tstack','curlF_snap','X_full','Y_full',"head")




%% get data stats
% fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',14);
% load(fname,'curlF')
% 
% dim_data= size(curlF);
% for xind = 1:dim_data(1)
%     for tind = 1:dim_data(3)
%         cFbr_data_var(xind,tind) = sqrt(var(curlF(xind,:,tind)));
%     end 
%     cFbr_data_stats(xind) = mean(cFbr_data_var(xind,:));
% end 
% head = 'data generated in generate_forcing_FFBL_v2 in the get data stats section';
%% %%%%%%%%%%%%%%%%%%% visual of the width function 

figure(10)
for tind = 1:nt
    subplot(121)
    plot(x,W_1s(1,:,tind))
    ylim([-15,15])
    ylabel('W(x-ct)*b(t)')
    xlabel('x(m)')
    cmocean('balance');
    caxis([-5,5])
    subplot(122)
    pcolorcen(X,Y,cFbr(:,:,tind))
    cmocean('balance');
    caxis([-0.3,0.3])
    colorbar
    title('$\nabla \times F_{br}$','Interpreter','latex')
    drawnow
    pause(0.2)
end 



%%
figure()
pcolorcen(-X,Y,cFbr(:,:,10))
cmocean('balance');
caxis([-0.3,0.3])


%% function 
function S_wb=make_WB(ky,ky0)
% make Weibull spectra using peak wavenumber ky0 and wavenumber vector ky
beta= 1.4; %modified gamma in the new paper version 
lam = ky0/((beta-1)/beta)^(1/beta); 
S_wb= (beta/(lam^beta))*abs(ky).^(beta-1).* exp(-(abs(ky/lam)).^(beta)); % Weibull spectra
end 