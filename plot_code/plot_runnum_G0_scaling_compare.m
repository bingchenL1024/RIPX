% Bingchen Liu June 14, 2024
% This code compare different scaling for each run
% Note this is in function format, and is included in master_plot_withloop
% to get plot for each/all run


function plot_runnum_G0_scaling_compare(ind)

    close all
    addpath('/data1/nkumar/RIPX/M_Files')
    addpath('/data1/bliu')
    
    load('/data1/bliu/data/RIPX_bath_guide.mat')
    load('/data1/bliu/data/SS_raw.mat')

    %% old exploratory code of using the full expression of derivative of Dw/c
    %plot_runnum_G0_scaling_full_derivative
    
    %% %%%%%%%%%%%%%%%%%%%%% collect data 
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

    %%%%%%%% Use Falk's code to get income wavenumber
    om = S(ind).fp_T.*(2.*pi);
    kw = real(get_wavenum(om,h2)); % use omega at breaking pt 
    kw_pm = real(kw./(2*pi));



%%  Calculate scaling parameters 
    
%%%%%%%%%%%%%%%%% d/dx(D_w)
%     dx=5;
%     [~,dDw]= calculate_gradient_from_central_difference(dECg,3,dx); %d/dx(D_w)
%     dDw = dDw';
%     dDw(end+1)= 0;
%     dDw= [0; dDw];
    
%%%%%%%%%%%%%%%% phase velocity calculation
    c_linear = real((g.* h2).^(1/2));
    c_nonlin= real((g.* h2).^(1/2).*((1-(kw.*h2).^2)./3).^(1/2));
    c_solit = real(sqrt(g.*(h2+Hs)));
 
%%%%%%%%%%%%%%%% Scaling using h (simplified)  
    %G0_scale_h_simp =(dECg./(real(sqrt(g.*h2.^3)).*h2)).*beta;
    G0_scale_h_simp =(dECg./(c_linear.*h2.^2));
    
%%%%%%%%%%%%%%% Scaling using kw
    G0_scale_k_h =(dECg./(c_linear.*h2)).*kw_pm;
    G0_scale_k_hb =(dECg./(c_linear).*hb).*kw_pm;

%%%%%%%%%%%%%%%%%%%%%% scaling using hb (at breaking) instead of h(x)
    G0_scale_hb_linear = (dECg./(c_linear.*h2.*hb));
    G0_scale_hb_nonlin = (dECg./(c_nonlin.*h2.*hb));
    G0_scale_hb_4all = (dECg./(c_linear.*hb.^2));
    
%%%%%%%%%%%%%%%%%%%% scaling using hb, h, and solitary wave phase velocity
    G0_scale_hb_solit = dECg./(c_solit.*h2.*hb);

%%%%%%%%%%%%%%%%%%%%%%%%% perturbation parameter calculation

    mu2= (kw.*h2).^2;
    epsl = Hs./h2;
    
%% %%%%%%%%%%%%%%%%   Plotting range setup
    xrange= 1:length(x2);
    yrange = 1:length(x2)-5;
    xdomain_comp = 1:length(x2);
    %xdomain_comp =length(x2)-10:length(x2);


    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% W/ hb
    %xrange=length(x2)-10:length(x2); %plot only section of the scaling to prevent sqrt(gh) blowing up 
    figure()
    sgtitle([runpar.wave,runpar.bath, ' w/ hb'],'Fontsize',20)

    subplot(3,2,1)
    plot(x,curlFstd,'Linewidth',3)
    hold on 
    scatter(x(ind_3loc_1),curlFstd(ind_3loc_1),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$std(\nabla \times F_{br}) \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(22)
    grid on   
    
    subplot(3,2,3)
    plot(x2(xrange),dECg(xrange),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),dECg(ind_3loc_2),80,"filled")
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac{d}{dx} \, E c_{g}$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(22)
    grid on     
    hold off 
    
    subplot(3,2,5)
    plot(x2(xrange),G0_scale_hb_4all(xrange),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_hb_4all(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$D_w/(c_l \, h_{b} ^2) \, (s^{-2})$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    %ylim([min(G0_scale_hb_4all(yrange)),max(G0_scale_hb_4all(yrange))])
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(22)
    grid on 
    
    subplot(3,2,2)
    plot(x2(xrange),mu2(xrange),'Linewidth',3)
    hold on 
    plot(x2(xrange),epsl(xrange),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),mu2(ind_3loc_2),80,"filled")
    hold on 
    scatter(x2(ind_3loc_2),epsl(ind_3loc_2),80,"filled")
    legend('$\mu ^2$','$\epsilon$','Interpreter','latex','Location','northwest')
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\mu^2 \hspace{0.1in} \mathrm{and} \hspace{0.1in} \epsilon$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    xlim([min(x2),max(x2)])
    ylim([min([mu2;epsl]),max([mu2;epsl])])
    niceplot_nobold_nomintick(22)
    grid on     
    hold off
    
    subplot(3,2,4)
    plot(x2(xrange),G0_scale_hb_linear(xrange),'Linewidth',3)
    hold on 
    plot(x2(xrange),G0_scale_hb_solit(xrange),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_hb_solit(ind_3loc_2),80,"filled")
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_hb_linear(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('${D_w}/({c \, h \, h_{b}}) \, (s^{-2})$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    legend('Linear c','Solitary c','Location','northwest')
    xlim([min(x2),max(x2)])
    ylim([min([G0_scale_hb_linear(yrange);G0_scale_hb_solit(yrange)]),max([G0_scale_hb_linear(yrange);G0_scale_hb_solit(yrange)])])
    niceplot_nobold_nomintick(22)
    grid on     
    
    
    subplot(3,2,6)
    plot(x2(xrange),G0_scale_k_hb(xrange),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_k_hb(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac {D_w}{c_l \, h_{b}} \, k_w \, (s^{-2})$','Interpreter','latex');
    %set(h,'Interpreter','latex')
    ylim([min(G0_scale_k_hb(yrange)),max(G0_scale_k_hb(yrange))])
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(22)
    grid on 


    width= 30;
    height = 15;
    set(gcf,'Units','inches','Position',[0,0,width,height])
    set(gcf,'visible','off')
    
    saveas(gcf,['/data1/bliu/figures/RIPX_allrun/scaling/with_hb/','run_',num2str(ind),'_scaling_w_hb.png'])
    
    
  
    
    
%% %%%%%%%%%%%%%%%%% G0 scaling using h (full) and with different terms     
%     xrange= 1:length(x2);
%     %xrange=length(x2)-10:length(x2); %plot only section of the scaling to prevent sqrt(gh) blowing up 
%     figure()
%     sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
% 
%     subplot(2,2,1)
%     plot(x,curlFstd,'Linewidth',3)
%     hold on 
%     scatter(x(ind_3loc_1),curlFstd(ind_3loc_1),80,"filled")
%     hold off 
%     xlabel('Cross-shore[x(m)]','Fontsize',16)
%     ylabel('$std(\nabla \times F_{br}) \, (s^{-2})$','Interpreter','latex')
%     xlim([min(x2),max(x2)])
%     niceplot_nobold_nomintick(20)
%     grid on   
%     
%     subplot(2,2,3)
%     plot(x2(xrange),G0_scale_h_full(xrange),'Linewidth',3)
%     hold on 
%     scatter(x2(ind_3loc_2),G0_scale_h_full(ind_3loc_2),80,"filled")
%     hold off 
%     xlabel('Cross-shore[x(m)]','Fontsize',16)
%     ylabel('$\frac{\frac{\partial}{\partial x} \, D_{w}}{c_{g} \, h} -\frac{3}{2} \frac{D_{w}}{c_{g} \, h^2} \, \beta \, (s^{-2})$','Interpreter','latex');
%     %set(h,'Interpreter','latex')
%     xlim([min(x2),max(x2)])
%     niceplot_nobold_nomintick(20)
%     grid on     
%     
%     subplot(2,2,2)
%     plot(x2(xrange),G0_1_scale(xrange),'Linewidth',3)
%     hold on 
%     scatter(x2(ind_3loc_2),G0_1_scale(ind_3loc_2),80,"filled")
%     hold off 
%     xlabel('Cross-shore[x(m)]','Fontsize',16)
%     ylabel('$\frac{\frac{\partial}{\partial x} \, D_{w}}{c_{g} \, h} \, (s^{-2})$','Interpreter','latex');
%     %set(h,'Interpreter','latex')
%     sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
%     xlim([min(x2),max(x2)])
%     niceplot_nobold_nomintick(20)
%     grid on     
%     
%     subplot(2,2,4)
%     plot(x2(xrange),G0_2_scale(xrange),'Linewidth',3)
%     hold on 
%     scatter(x2(ind_3loc_2),G0_2_scale(ind_3loc_2),80,"filled")
%     hold off 
%     xlabel('Cross-shore[x(m)]','Fontsize',16)
%     ylabel('$\frac{3}{2} \frac{D_{w}}{c_{g} \, h^2} \, \beta \, (s^{-2})$','Interpreter','latex');
%     %set(h,'Interpreter','latex')
%     sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
%     xlim([min(x2),max(x2)])
%     niceplot_nobold_nomintick(20)
%     grid on 

    





end 
