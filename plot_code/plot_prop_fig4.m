% Bingchen Liu April, 2024
% This code is used to plot proposal Fig 4 (c) and (d) 
% Note it requires user input of specific wavefield and bath


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


%% Specific wavefield (requires user input)

dir_ds = '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
varname_ds = 'runlist.amp.amp_';
dir_bath= '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
varname_bath= 'runlist.bath.bath_';


% index with specific wave height and peak period, NOT general 
    input_Hs= sprintf('%.1f',input('Enter the significant wave height you want (e.g.1.1): '));
    input_Tp= sprintf('%.1f',input('Enter the peak period you want (e.g. 14.0): '));
    input_bath= input('Enter the slope of the bathymetry you want (e.g. 004):','s');

    ind_2_5_spec=data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds2.5'])';
    ind_5_spec=data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds5.0'])';
    ind_10_spec=data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds10.0'])';
    ind_20_spec=data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds20.0'])';
    
    ind_bath_temp=data_filter(dir_bath,varname_bath,input_bath)';
    
    ind_2_5_bath_ds_specific = ind_2_5_spec(ismember(ind_2_5_spec,ind_bath_temp));
    ind_5_bath_ds_specific = ind_5_spec(ismember(ind_5_spec,ind_bath_temp));
    ind_10_bath_ds_specific = ind_10_spec(ismember(ind_10_spec,ind_bath_temp));
    ind_20_bath_ds_specific = ind_20_spec(ismember(ind_20_spec,ind_bath_temp));
    ind = [ind_2_5_bath_ds_specific,ind_5_bath_ds_specific,ind_10_bath_ds_specific,ind_20_bath_ds_specific]; %index for specific plot use, include specific bathy


%% special case 
for i = 1:4
    filename1 = ['/data1/nkumar/RIPX/M_Files/RIPX_mstd_vort/RIPX_mstd_vort_',sprintf('%04d',ind(i)),'.mat'];
    load(filename1)
    std_vort_ind{i} = vort_mstd;
    clear vort_mstd
end 

for i = 1:4
    x_ind{i} = x{ind(i)};   
    x2_ind{i} = x2{ind(i)};
    h_ind{i} = h{ind(i)};    
    h2_ind{i} = h2{ind(i)};
    std_curlF_ind{i} = curlFstd{ind(i)};
end

%% plot 
figure()
subplot(3,1,1)
plot(x_ind{1},h_ind{1},'Linewidth',3)
hold on 
plot(x_ind{2},h_ind{2},'Linewidth',3)
hold on 
plot(x_ind{3},h_ind{3},'Linewidth',3)
hold on 
plot(x_ind{4},h_ind{4},'Linewidth',3)
legend('\sigma_\theta = 2.5^{\circ}','\sigma_\theta = 5^{\circ}','\sigma_\theta = 10^{\circ}','\sigma_\theta = 20^{\circ}','location','northwest')
xlabel('Cross-shore[x(m)]','Fontsize',16)
ylabel('Water depth (m)','Fontsize',16)
set(gca,'ydir','reverse')
niceplot_nobold_nomintick(18)
hold off
grid on 

subplot(3,1,2)
plot(x_ind{1},std_vort_ind{1},'Linewidth',3)
hold on 
plot(x_ind{2},std_vort_ind{2},'Linewidth',3)
hold on 
plot(x_ind{3},std_vort_ind{3},'Linewidth',3)
hold on 
plot(x_ind{4},std_vort_ind{4},'Linewidth',3)
legend('\sigma_\theta = 2.5^{\circ}','\sigma_\theta = 5^{\circ}','\sigma_\theta = 10^{\circ}','\sigma_\theta = 20^{\circ}','location','northwest')
xlabel('Cross-shore[x(m)]','Fontsize',16)
ylabel({'Vorticity','std(\omega)(s^{-1})'},'Fontsize',16)
niceplot_nobold_nomintick(18)
hold off
grid on

subplot(3,1,3)
plot(x_ind{1},std_curlF_ind{1},'Linewidth',3)
hold on 
plot(x_ind{2},std_curlF_ind{2},'Linewidth',3)
hold on 
plot(x_ind{3},std_curlF_ind{3},'Linewidth',3)
hold on 
plot(x_ind{4},std_curlF_ind{4},'Linewidth',3)
legend('\sigma_\theta = 2.5^{\circ}','\sigma_\theta = 5^{\circ}','\sigma_\theta = 10^{\circ}','\sigma_\theta = 20^{\circ}','location','northwest')
xlabel('Cross-shore[x(m)]','Fontsize',16)
ylabel({'Vorticity Forcing','std(\nabla \times F_{br})(s^{-2})'},'Fontsize',16)
niceplot_nobold_nomintick(18)
%sgtitle(['Hs = ',input_Hs,' m',', Tp = ',input_Tp,' s',', beach slope = ', input_bath],'Fontsize',20)
sgtitle(['Hs = ',input_Hs,' m',', Tp = ',input_Tp,' s',', beach slope = 0.04' ],'Fontsize',20)
grid on 
hold off