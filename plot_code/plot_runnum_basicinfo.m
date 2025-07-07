%Bingchen Liu April 2024
%This code is used to obtain plot after specifying RUN NUMBER INDEX INSTEAD of
%run parameter (use plot_onerun for that function)
% Note it only plot one run

function  plot_runnum_basicinfo(ind)
    close all
    addpath('/data1/nkumar/RIPX/M_Files')
    addpath('/data1/bliu')
    
    load('/data1/bliu/data/RIPX_bath_guide.mat')
    load('/data1/bliu/data/SS_raw.mat')

%     dir_ds = '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
%     varname_ds = 'runlist.amp.amp_';
%     dir_bath= '/data1/nkumar/RIPX/M_Files/Important_Files/RIPX_bath_guide.mat';
%     varname_bath= 'runlist.bath.bath_';
% 
%     input_Hs= sprintf('%.1f',input('Enter the significant wave height you want (e.g.1.1): '));
%     input_Tp= sprintf('%.1f',input('Enter the peak period you want (e.g. 14.0): '));
%     input_dr= sprintf('%.1f',input('Enter the directional spread you want (e.g. 2.5): '));
%     input_bath= input('Enter the slope of the bathymetry you want (e.g. 002-006):','s');
%     
%     ind_wave =data_filter(dir_ds,varname_ds,['Hs',input_Hs,'_Tp',input_Tp,'_ds',input_dr])';
%     ind_bath2=data_filter(dir_bath,varname_bath,input_bath)';
%     
%     ind = ind_wave(ismember(ind_wave,ind_bath2));
    
%% %%%%%%%%%%%%%%%%%%%%% collect data 
    ind_3loc_1 = get_3locs(S(ind),1);
    ind_3loc_2 = get_3locs(S(ind),2);
    
    runpar= get_runpara(ind);

    g=9.81;

    
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
    
    filename1 = ['/data1/nkumar/RIPX/M_Files/RIPX_mstd_vort/RIPX_mstd_vort_',sprintf('%04d',ind),'.mat'];
    load(filename1)
    std_vort = vort_mstd;   

    beta= AA.A(ind,1); %beach slope

    %%%%%%%% Use Falk's code to get income wavenumber
    om = S(ind).fp_T.*(2.*pi);
    kw = get_wavenum(om,h2); % use omega at breaking pt 
    kw_pm = (kw./(2*pi));

    %%%%%%%%%%%%%%%%%%%%% calculate d/dx(D_w)
    dx=5;
    [~,dDw]= calculate_gradient_from_central_difference(dECg,3,dx); %d/dx(D_w)
    dDw = dDw';
    dDw(end+1)= 0;
    dDw= [0; dDw];


%%  Plotting    
%%%%%%%%%%%%%%%%%%%%%% profile of general quantity
    figure()
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20)

    subplot(3,3,1)
    plot(x2,h2,'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),h2(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('Water depth (m)','Fontsize',16)
    set(gca,'ydir','reverse')
    xlim([min(x2),max(x2)])
    grid on 
    niceplot_nobold_nomintick(20)
    
    
    subplot(3,3,4)
    plot(x2,Hs,'Linewidth',3)
    hold on
    scatter(x2(ind_3loc_2),Hs(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('Hs (m)','Fontsize',16)
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 
    
    subplot(3,3,7)
    plot(x2,sigma_th,'Linewidth',3)
    hold on
    scatter(x2(ind_3loc_2),sigma_th(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('\sigma_{\theta}','Fontsize',16)
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    subplot(3,3,2)
    plot(x2,ECg,'Linewidth',3)
    hold on
    scatter(x2(ind_3loc_2),ECg(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('ECg (J m^{-1} s^{-1})','Fontsize',16)
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 

    subplot(3,3,5)
    plot(x2,dECg,'Linewidth',3)
    hold on
    scatter(x2(ind_3loc_2),dECg(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac{d}{dx} ECg (m^{3} s^{-3})$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 
    
    subplot(3,3,8)
    plot(x2,dDw,'Linewidth',3)
    hold on
    scatter(x2(ind_3loc_2),dDw(ind_3loc_2),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac{d}{dx} D_w (m^{2} s^{-3})$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,3)
    plot(x,Fbr_x,'Linewidth',3)
    hold on 
    plot(x,Fbr_y,'Linewidth',3)
    hold on 
    plot(x,Fbr_mag,'Linewidth',3)
    hold on
    scatter(x(ind_3loc_1),Fbr_mag(ind_3loc_1),80,"filled")
    legend('x-component','y-component','Norm','3 sz loc','Location','northwest')
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$F_{br} (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on
    hold off 
    
    subplot(3,3,6)
    plot(x,std_vort,'Linewidth',3)
    hold on
    scatter(x(ind_3loc_1),std_vort(ind_3loc_1),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('std(\omega)(s^{-1})')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 

    subplot(3,3,9)
    plot(x,curlFstd,'Linewidth',3)
    hold on
    scatter(x(ind_3loc_1),curlFstd(ind_3loc_1),80,"filled")
    hold off
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('std(\nabla \times F_{br})(s^{-2})')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 

    width= 30;
    height = 15;
    set(gcf,'Units','inches','Position',[0,0,width,height])
    set(gcf,'visible','off') 
  
    saveas(gcf,['/data1/bliu/figures/RIPX_allrun/basic_info/version2/','run_',num2str(ind),'.png'])
end 