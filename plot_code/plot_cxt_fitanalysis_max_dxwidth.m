% Bingchen Liu, Mar 3, 2025
% This code analyze the fit parameter a that from data without
% interpolation
% only include the good ones for paper 
clear
load('/data1/bliu/data/cxt_fitanalysis_max_dxwidth.mat')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_nointerp_runinfo_G0lim0p15.mat')

load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');

runnum = load('runnum_slp.mat');


%% July 9 new plot 
% plot individual run to see overall pattern 
close all hidden
figure()
%subplot(211)
for runind = 1:24
        clf
        subplot(311)
        runind
        scatter(fitpara.slp2(runind).x_nond,fitpara.slp2(runind).a,60,'r','filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp3(runind).x_nond,fitpara.slp3(runind).a,60,'g','filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp4(runind).x_nond,fitpara.slp4(runind).a,60,'b','filled','Marker','^','MarkerEdgeColor','k')
        hold off
        ylim([0.3,2])
        runparaslp2 = get_runpara(runnum.runnum_slp.slp2(runind))
        runparaslp3 = get_runpara(runnum.runnum_slp.slp3(runind))
        runparaslp4 = get_runpara(runnum.runnum_slp.slp4(runind))
        % hb_slp2 = S(runnum.runnum_slp.slp2(runind)).hb
        % hb_slp3 = S(runnum.runnum_slp.slp3(runind)).hb
        % hb_slp4 = S(runnum.runnum_slp.slp4(runind)).hb
        xlabel('x/L_{sz}','FontSize',26)
        ylabel('$\tau$ (s)','interpreter','latex','FontSize',26)
        title(runparaslp2.wave,'FontSize',30)
        subplot(312)
        scatter(fitpara.slp2(runind).x_nond,fitpara.slp2(runind).G0,60,'r','filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp3(runind).x_nond,fitpara.slp3(runind).G0,60,'g','filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp4(runind).x_nond,fitpara.slp4(runind).G0,60,'b','filled','Marker','^','MarkerEdgeColor','k')
        hold off
        xlabel('x/L_{sz}','FontSize',26)
        ylabel('$G_{0} (m^{-2})$','interpreter','latex','FontSize',26)

        subplot(313)
        scatter(fitpara.slp2(runind).x_nond,fitpara.slp2(runind).Dw,60,'r','filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp3(runind).x_nond,fitpara.slp3(runind).Dw,60,'g','filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp4(runind).x_nond,fitpara.slp4(runind).Dw,60,'b','filled','Marker','^','MarkerEdgeColor','k')
        hold off
        xlabel('x/L_{sz}','FontSize',26)
        ylabel('$D_{w} (m^3 s^{-3})$','interpreter','latex','FontSize',26)

        width= 15;
        height = 15;
        set(gcf,'Units','inches','Position',[0,0,width,height])
        set(gcf,'visible','off') 
  
        saveas(gcf,['/data1/bliu/figures/RIPX_allrun/tau/','run_',num2str(runind),'.png'])
end 


%% July 9 plot tau vs skill 
% plot individual run to see overall pattern 
figure()
%subplot(211)
for runind = 1:24
        
        scatter(fitpara.slp2(runind).rsq,fitpara.slp2(runind).a,60,'r','filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp3(runind).rsq,fitpara.slp3(runind).a,60,'g','filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara.slp4(runind).rsq,fitpara.slp4(runind).a,60,'b','filled','Marker','^','MarkerEdgeColor','k')
        hold on
        
end 
hold off 
xlabel('fit skill','FontSize',26)
ylabel('$\tau$ (s)','Interpreter','latex','FontSize',26)
%% tau/\sqrt(h/g) VS nond G0 
col = cmocean('thermal',5);

figure()
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).G0_nond(xind),fitpara_5loc.slp2(runind).a(xind)/fitpara_5loc.slp2(runind).t_scale(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp3(runind).G0_nond(xind),fitpara_5loc.slp3(runind).a(xind)/fitpara_5loc.slp3(runind).t_scale(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp4(runind).G0_nond(xind),fitpara_5loc.slp4(runind).a(xind)/fitpara_5loc.slp4(runind).t_scale(xind),50,col(xind,:),'filled')
        hold on 
    end 
end 
%legend('slp2','slp3','slp4')
xlabel('${G_0}/{G_{\mathrm{max}}}$','Interpreter','latex')
ylabel('$\tau / \sqrt{h/g}$','Interpreter','latex')
niceplot(20)
hold off 

%% nondG0 VS x/xb (colored by slp) 
col = cmocean('thermal',5);

figure()
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).x_nond(xind),fitpara_5loc.slp2(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp3(runind).x_nond(xind),fitpara_5loc.slp3(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp4(runind).x_nond(xind),fitpara_5loc.slp4(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
    end 
end 
%legend('slp2','slp3','slp4')
ylabel('${G_0}/{G_{\mathrm{max}}}$','Interpreter','latex')
xlabel('$x/x_{\mathrm{b}}$','Interpreter','latex')
niceplot(20)
hold off 

%% nondG0 VS h (colored by slp) 
col = cmocean('thermal',5);

figure()
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).h(xind),fitpara_5loc.slp2(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp3(runind).h(xind),fitpara_5loc.slp3(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
        scatter(fitpara_5loc.slp4(runind).h(xind),fitpara_5loc.slp4(runind).G0_nond(xind),50,col(xind,:),'filled')
        hold on 
    end 
end 
%legend('slp2','slp3','slp4')
ylabel('${G_0}/{G_{\mathrm{max}}}$','Interpreter','latex')
xlabel('$h(x)$','Interpreter','latex')
niceplot(20)
hold off 
%% tau(s)  histogram

figure()
subplot(121)
histogram(a_tot_all,'Normalization','probability')
xlabel('$\tau$','Interpreter','latex')
ylabel('$PDF=c_{i}/N_{\mathrm{tot}}$','Interpreter','latex')
title('Every 1 m')
niceplot(22)

subplot(122)
histogram(t_tot_5loc,'Normalization','probability')
xlabel('$\tau$','Interpreter','latex')
ylabel('$PDF=c_{i}/N_{\mathrm{tot}}$','Interpreter','latex')
title('5 locations in crosshore')
niceplot(22)
%% tau(s) VS x/xb (colored by slp) 
figure()
% scatter(t_scale_all,a_tot_all,45,'filled')
hold on 
scatter(x_nond_tot.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(x_nond_tot.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(x_nond_tot.slp4,a_tot.slp4,45,'filled')
hold off
xlabel('$\frac{x}{x_{b}}$','Interpreter','latex')
ylabel('$\tau (s)$','Interpreter','latex')
niceplot(18)

%% tau(s) VS sqrt(h/g) (colored by slp) 
figure()
% scatter(t_scale_all,a_tot_all,45,'filled')
hold on 
scatter(t_scale.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(t_scale.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(t_scale.slp4,a_tot.slp4,45,'filled')
hold off
xlabel('$\sqrt{h/g} (s)$','Interpreter','latex')
ylabel('$\tau (s)$','Interpreter','latex')
%title('$cxt_{fit} = \exp{(-\frac{x}{a})} \cos{(\frac{x}{b} +c)}/\cos{c}$','interpreter','latex')
%title({['a/h = ',num2str(f_a_nond.A),'*Hs/h +',num2str(f_a_nond.B)],[' r^2 = ',num2str(gof_a_nond.rsquare)]},'interpreter','latex')
%legend('slp2','slp3','slp4')
niceplot(18)
xlim([0,Inf])
ylim([0,Inf])

%% tau(s) VS directional spread (colored by slp) 
figure()
% scatter(t_scale_all,a_tot_all,45,'filled')
hold on 
scatter(dirspr.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(dirspr.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(dirspr.slp4,a_tot.slp4,45,'filled')
hold off
xlabel('$\sigma_{\theta b}$','Interpreter','latex')
ylabel('$\tau (s)$','Interpreter','latex')
niceplot(18)
%% nond t vs nond G0
figure()
scatter(G0_nond.slp2,a_tot.slp2./t_scale.slp2,50,'filled')
hold on 
scatter(G0_nond.slp3,a_tot.slp3./t_scale.slp3,50,'filled')
hold on 
scatter(G0_nond.slp4,a_tot.slp4./t_scale.slp4,50,'filled')
hold on 
%plot(nond_G0_model,a_nond_model_G0,'LineWidth',6)
ylim([0,10])
%title(['r^2=',num2str(gof_a_nond_G0.rsquare)])
legend('slp2','slp3','slp4')
xlabel('$\frac{G_0}{G_{0max}}$','Interpreter','latex')
ylabel('$\frac{\tau}{\sqrt{h/g}}$','Interpreter','latex')
niceplot(20)
hold off 


%% nond t vs dirspr
figure()
scatter(dirspr.slp2,a_tot.slp2./t_scale.slp2,50,'filled')
hold on 
scatter(dirspr.slp3,a_tot.slp3./t_scale.slp3,50,'filled')
hold on 
scatter(dirspr.slp4,a_tot.slp4./t_scale.slp4,50,'filled')
hold on 
%plot(nond_G0_model,a_nond_model_G0,'LineWidth',6)
ylim([0,10])
%title(['r^2=',num2str(gof_a_nond_G0.rsquare)])
legend('slp2','slp3','slp4')
xlabel('$\sigma_{\theta b}$','Interpreter','latex')
ylabel('$\frac{\tau}{\sqrt{h/g}}$','Interpreter','latex')
niceplot(20)
hold off 
%% t(s) vs 1/G0(s) -- 5 loc

figure()
for runind = 1:24
scatter(G0.G0_data.slp2(runind,:),fitpara_5loc.slp2(runind).a,50,'filled','k')
hold on 
scatter(G0.G0_data.slp3(runind,:),fitpara_5loc.slp3(runind).a,50,'filled','k')
hold on 
scatter(G0.G0_data.slp4(runind,:),fitpara_5loc.slp4(runind).a,50,'filled','k')
hold on 
end 
%legend('slp2','slp3','slp4')
xlabel('${G_0}$','Interpreter','latex')
ylabel('$\tau$','Interpreter','latex')
niceplot(20)
hold off 


