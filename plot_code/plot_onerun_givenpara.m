%Bingchen Liu March 2024
%This code is used to obtain plot after specifying run parameters (e.g. Hs, Tp)
% Note it only plot one run

function plot_onerun()
    addpath('/data1/nkumar/RIPX/M_Files')
    addpath('/data1/bliu')
    
    load('/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat')
    load('/data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat')

    dir_ds = '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
    varname_ds = 'runlist.amp.amp_';
    dir_bath= '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
    varname_bath= 'runlist.bath.bath_';
    
    input_Hs= sprintf('%.1f',input('Enter the significant wave height you want (e.g. 0.5; 0.8; 1.1): '));
    input_Tp= sprintf('%.1f',input('Enter the peak period you want (e.g. 8.0; 14.0; ): '));
    input_dr= sprintf('%.1f',input('Enter the directional spread you want (e.g. 2.5; 5.0; 10.0; 20.0): '));
    input_bath= input('Enter the slope of the bathymetry you want (e.g. 002-006):','s');
    
    ind_wave =data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds',input_dr])';
    ind_bath2=data_filter(dir_bath,varname_bath,input_bath)';
    
    ind = ind_wave(ismember(ind_wave,ind_bath2));
    
    x= S(ind).X;     
    x2 = S(ind).X2;
    h = S(ind).h;
    h2=S(ind).h2;
    curlFstd = S(ind).curlF_std;
    Hs = S(ind).Hs;
    sigma_th = S(ind).sigma_th;
    mean_th = S(ind).mean_th;
    
    
    filename1 = ['/data1/nkumar/RIPX/M_Files/RIPX_mstd_vort/RIPX_mstd_vort_',sprintf('%04d',ind),'.mat'];
    load(filename1)
    std_vort = vort_mstd;
    
    figure()
    subplot(3,2,1)
    plot(x2,h2,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('Water depth (m)','Fontsize',16)
    
    subplot(3,2,2)
    plot(x2,Hs,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('Significant wave height (m)','Fontsize',16)
    
    subplot(3,2,3)
    plot(x2,sigma_th,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('\sigma_{\theta}','Fontsize',16)
    
    
    subplot(3,2,4)
    plot(x2,mean_th,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel(' mean {\theta}','Fontsize',16)
    
    subplot(3,2,5)
    plot(x,std_vort,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('std(\omega)(s^{-1})','Fontsize',16)

    subplot(3,2,6)
    plot(x,curlFstd,'Linewidth',3)
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('std(\nabla \times F_{br})(s^{-2})','Fontsize',16)
    
    sgtitle(['Hs = ',input_Hs,', Tp = ',input_Tp, ', \sigma_{\theta} = ',input_dr,', bathy:',input_bath],'Fontsize',20)
end 