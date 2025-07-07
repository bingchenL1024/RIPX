
% Bingchen Liu April, 2024
% This code is used to plot along-shore wavenumber spectra 
% Correspond to Fig 5 without the Weibull fit,
% NOTE: this code uses specific run, for general spec plot, check 

clear
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
addpath(genpath('/data1/bliu'))

load('/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat')
load('/data1/nkumar/RIPX/M_Files/Important_Files/bathRIPX2.mat')
load('/data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat')


%% load all the variabes needed 

for i = 1:120
    x{i}= S(i).X;     
    x2{i} = S(i).X2;
    h{i} = S(i).h;
    h2{i}=S(i).h2;
    curlFstd{i} = S(i).curlF_std;
    hwater{i} = S(i).h; 
    hwater2{i} = S(i).h2; 
    H_sig{i} = S(i).Hs;
    Ecg{i} = -S(i).ECG;
    dEcg{i} = -S(i).dECG;
    x_br{i} = S(i).xb; %Note that x_br is always positive. If want to use it to find ind_br, use NEGATIVE x_br
    Dw{i}=abs(dEcg{i})';%.*1020.*2;
    ds{i}= S(i).sigma_th;
    %[~,ind_xb{i}] = min(abs(-x_br{i}-x2{i}));
    ind_xb{i}= find(x2{i}==-x_br{i}); %find the index in x2 of breaking point 
    fp_br{i} = S(i).fp_T;
end 


%% constant input
g=9.81; %gravitational constant 
sz_crossloc= [0.33,0.5,0.75];%index for surf zone cross-shore location 



%% filter run with given parameters 



dir_ds = '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
varname_ds = 'runlist.amp.amp_';
dir_bath= '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
varname_bath= 'runlist.bath.bath_';

ind_run = 1:120;


%General index, ONLY include ds, OR bath, OR Hs, OR Tp
ind_bath2=data_filter(dir_bath,varname_bath,'002')';
ind_bath3=data_filter(dir_bath,varname_bath,'003')';
ind_bath4=data_filter(dir_bath,varname_bath,'004')';
ind_bath5=data_filter(dir_bath,varname_bath,'005')';
ind_bath6=data_filter(dir_bath,varname_bath,'006')';
ind_bath = [ind_bath2,ind_bath3,ind_bath4,ind_bath5,ind_bath6];%index for bathy for general use

ind_2_5_general=data_filter(dir_ds,varname_ds,'ds2.5')';
ind_5_general=data_filter(dir_ds,varname_ds,'ds5.0')';
ind_10_general=data_filter(dir_ds,varname_ds,'ds10.0')';
ind_20_general=data_filter(dir_ds,varname_ds,'ds20.0')';
ind_ds = [ind_2_5_general,ind_5_general,ind_10_general,ind_20_general]; %index for ds for general use 


spec = load('/data1/nkumar/RIPX/M_Files/RIPX_mean_curlF_kf_mspec/RIPX_mean_curlF_kf_mspec_0043.mat');
ind_spectra_all=data_filter(dir_ds,varname_ds,'Hs0.8_Tp8.0_ds10.0');
ind_spectra= ind_spectra_all((ismember(ind_spectra_all,ind_bath4)));

[h_max,ind_hmax] = max(H_sig{ind_spectra});
x_sz_exact = x2{ind_spectra}(ind_hmax);
[x_75sz,ind_75sz] = min(abs(x{ind_spectra}-x_sz_exact*sz_crossloc(3)));
[x_50sz,ind_50sz] = min(abs(x{ind_spectra}-x_sz_exact*sz_crossloc(2)));
[x_33sz,ind_33sz] = min(abs(x{ind_spectra}-x_sz_exact*sz_crossloc(1)));


figure()
loglog(spec.fck/(2*pi),spec.Sckm(ind_75sz,:),'Linewidth',3)
hold on 
loglog(spec.fck/(2*pi),spec.Sckm(ind_50sz,:),'Linewidth',3)
hold on 
loglog(spec.fck/(2*pi),spec.Sckm(ind_33sz,:),'Linewidth',3)
xlabel('ky(cpm)','FontSize',16)
ylabel('S_{\nabla \times F_{br}}(s^{-4}/cpm)','FontSize',16)
legend('x/L_{sz} = -0.75','x/L_{sz} = -0.50','x/L_{sz} = -0.33','FontSize',13)
axis([1e-3,0.2 1e-4 1e-2])
title('Hs = 0.8 m; Tp = 8.0 s','FontSize',20)
grid on 
hold off   