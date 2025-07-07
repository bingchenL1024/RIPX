%Bingchen Liu, July 10, 2024
%old version when exploring different scaling of G0 using full derivative
%expression
%Note: it might missing the function
%'calculate_gradient_from_central_difference'



function plot_runnum_G0_scaling_full_derivative(ind)

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
%%%%%%%%%%%%%%%% Fbr scaling    
    Fbr_scale= dECg./(real(sqrt(g.*h2)).*h2); %Fbr scaling from LH (1970) and Smith(2006), depth avg
    Fbr_scale= dECg./(c_linear.*hb);
%%%%%%%%%%%%%%% Scaling using h (full)   
    G0_2_scale = (3/2).*beta.*(dECg./(real(sqrt(g.*h2.^3)).*h2));%one term of curlFbr from chain rule 
    G0_1_scale = dDw./(real(sqrt(g.*h2.^3)));
    
    G0_2_scale = (3/2).*beta.*(dECg./(c_linear.*h2.^2));
    G0_1_scale = dDw./(c_linear.*h2);%another term of curlFbr from chain rule 
    G0_scale_h_full= G0_1_scale - G0_2_scale;
    
    
%%%%%%% %%%%%%%%%%%%%%%%%%% W/o hb
    % G0 scaling comparison with h instead of hb(detailed terms from chain rule -- maybe not very useful)
    close all
    figure()
    sgtitle([runpar.wave,runpar.bath, ' w/o hb'],'Fontsize',20)

    
    subplot(4,1,1)
    plot(x,curlFstd,'Linewidth',3)
    hold on 
    scatter(x(ind_3loc_1),curlFstd(ind_3loc_1),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$std(\nabla \times F_{br}) \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on
    
    
    subplot(4,1,2)
    plot(x2(xdomain_comp),G0_scale_h_simp(xdomain_comp),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_h_simp(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac{D_{w}}{c_{g} \, h} \, \frac{1}{h} \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on  
    
    subplot(4,1,3)
    plot(x2(xdomain_comp),G0_scale_h_full(xdomain_comp),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_h_full(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$\frac{\frac{\partial}{\partial x} \, D_{w}}{c_{g} \, h} -\frac{3}{2} \frac{D_{w}}{c_{g} \, h^2} \, \beta \, (s^{-2})$','Interpreter','latex');
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on 
    
    
    subplot(4,1,4)
    plot(x2(xdomain_comp),G0_scale_k_h(xdomain_comp),'Linewidth',3)
    hold on 
    scatter(x2(ind_3loc_2),G0_scale_k_h(ind_3loc_2),80,"filled")
    hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$ \frac{D_{w}}{c_{g} \, h} \, k_w \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(20)
    grid on      
 

    width= 30;
    height = 15;
    set(gcf,'Units','inches','Position',[0,0,width,height])
    set(gcf,'visible','off')
    
    %saveas(gcfc,['/data1/bliu/figures/RIPX_allrun/scaling/without_hb/','run_',num2str(ind),'_scaling_wo_hb.png'])