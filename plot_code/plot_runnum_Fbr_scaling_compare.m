



function plot_runnum_Fbr_scaling_compare(ind)

    close all
    addpath('/data1/nkumar/RIPX/M_Files')
    addpath(genpath('/data1/bliu'))
    
    load('/data1/bliu/data/RIPX_bath_guide.mat')
    load('/data1/bliu/data/SS_raw.mat')

    
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
    kw = get_wavenum(om,h2); % use omega at breaking pt 
    kw_pm = (kw./(2*pi));



%%  Calculate scaling parameters 
    
%%%%%%%%%%%%%%%%% d/dx(D_w)
    dx=5;
    [~,dDw]= calculate_gradient_from_central_difference(dECg,3,dx); %d/dx(D_w)
    dDw = dDw';
    dDw(end+1)= 0;
    dDw= [0; dDw];
    
%%%%%%%%%%%%%%%% phase velocity calculation
    c_linear = real((g.* h2).^(1/2));
    c_nonlin= real((g.* h2).^(1/2).*((1-(kw.*h2).^2)./3).^(1/2));
%%%%%%%%%%%%%%%% Fbr scaling    
    Fbr_scale_h= dECg./(c_linear.*h2); %Fbr scaling from LH (1970) and Smith(2006), depth avg
    Fbr_scale_hb= dECg./(c_linear.*hb);
%%%%%%%%%%%%%%%%%%%%%%%%% perturbation parameter calculation

    mu2= (kw.*h2).^2;
    epsl = Hs./h2;
    
%% %%%%%%%%%%%%%%%%   Plotting range setup
    xrange= 1:length(x2);
    yrange = 1:length(x2)-10;
    xdomain_comp = 1:length(x2);
    %xdomain_comp =length(x2)-10:length(x2);
%% plot 

    figure()
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20)

    subplot(3,1,1)
    plot(x,Fbr_mag,'Linewidth',3)
    hold on 
    scatter(x(ind_3loc_1),Fbr_mag(ind_3loc_1),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$|F_{br}| \, (m \, s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on

    subplot(3,1,2)
    plot(x2(xdomain_comp),Fbr_scale_hb(xdomain_comp),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),Fbr_scale_hb(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$F_{br} \sim \frac{D_{w}}{c_{g} \, h_b} \, (m \, s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 
    
    subplot(3,1,3)
    plot(x2(xdomain_comp),Fbr_scale_h(xdomain_comp),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),Fbr_scale_h(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$F_{br} \sim \frac{D_{w}}{c_{g} \, h} \, (m \, s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 

    
    width= 15;
    height = 15;
    set(gcf,'Units','inches','Position',[0,0,width,height])
    set(gcf,'visible','off')
    
    saveas(gcf,['/data1/bliu/figures/RIPX_allrun/Fbr/','run_',num2str(ind),'.png'])

end 