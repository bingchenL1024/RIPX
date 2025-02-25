% Bingchen Liu June 22, 2024
% Pending! might not be very useful and will incorporate into plot code due
% to index input requirement June 22, 2024
% This is code calculate the profile of different scaling
% The data output from this code is used in "plot_runnum_scaling_compare.m"
% function



close all
addpath('/data1/nkumar/RIPX/M_Files')
addpath('/data1/bliu')

load('/data1/bliu/data/RIPX_bath_guide.mat')
load('/data1/bliu/data/SS_raw.mat')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% collect data 
g=9.81;
    
ind_3loc_1 = get_3locs(S(ind),1);
ind_3loc_2 = get_3locs(S(ind),2);

    
runpar= get_runpara(ind);

x= S(ind).X;     
x2 = S(ind).X2;
h = S(ind).h;
h2=S(ind).h2;
curlFstd = S(ind).curlF_std;
Hs = S(ind).Hs;
sigma_th = S(ind).sigma_th;
mean_th = S(ind).mean_th;
ECg= -real(S(ind).ECG);
dECg =real(S(ind).dECG');
Fbr_mag = S(ind).Fbr_mag;
Fbr_x = S(ind).Fbr_x;
Fbr_y = S(ind).Fbr_y;
hb= S(ind).hb;
    
filename1 = ['/data1/nkumar/RIPX/M_Files/RIPX_mstd_vort/RIPX_mstd_vort_',sprintf('%04d',ind),'.mat'];
load(filename1)
std_vort = vort_mstd;   

beta= AA.A(ind,1); %beach slope

%% Use Falk's code to get income wavenumber
om = S(ind).fp_T.*(2.*pi);
kw = get_wavenum(om,h2); % use omega at breaking pt 
kw_pm = (kw./(2*pi));



Calculate scaling parameters 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% d/dx(D_w)
dx=5;
[~,dDw]= calculate_gradient_from_central_difference(dECg,3,dx); %d/dx(D_w)
dDw = dDw';
dDw(end+1)= 0;
dDw= [0; dDw];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% phase velocity calculation
c_linear = real((g.* h2).^(1/2));
c_nonlin= real((g.* h2).^(1/2).*((1-(kw.*h2).^2)./3).^(1/2));
%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%% Fbr scaling    
Fbr_scale= dECg./(real(sqrt(g.*h2)).*h2); %Fbr scaling from LH (1970) and Smith(2006), depth avg
Fbr_scale= dECg./(c_linear.*h2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling using h (full)   
G0_2_scale = (3/2).*beta.*(dECg./(real(sqrt(g.*h2.^3)).*h2));%one term of curlFbr from chain rule 
G0_1_scale = dDw./(real(sqrt(g.*h2.^3)));
    
G0_2_scale = (3/2).*beta.*(dECg./(c_linear.*h2.^2));
G0_1_scale = dDw./(c_linear.*h2);%another term of curlFbr from chain rule 
G0_scale_h_full= G0_1_scale - G0_2_scale;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling using h (simplified)  
G0_scale_h_simp =(dECg./(real(sqrt(g.*h2.^3)).*h2)).*beta;
G0_scale_h_simp =(dECg./(c_linear.*h2.^2)).*beta;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scaling using kw
G0_scale_k =(dECg./(real(sqrt(g.*h2.^3)))).*kw_pm.*beta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scaling using hb (at breaking) instead of h(x)
G0_scale_hb_linear = (dECg./(c_linear.*h2.*hb)).*beta;
G0_scale_hb_nonlin = (dECg./(c_nonlin.*h2.*hb)).*beta;
    
G0_scale_hb_4all = (dECg./(c_linear.*hb.^2)).*beta;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% perturbation parameter calculation

mu2= (kw.*h2).^2;
epsl = Hs./h2;
    