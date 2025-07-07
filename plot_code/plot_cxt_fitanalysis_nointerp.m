% Bingchen Liu, Nov 14, 2024
% This code analyze the fit parameter a that from data without
% interpolation
% only include the good ones for paper 
clear
load('/data1/bliu/data/cxt_nointerp_fitanalysis')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_nointerp_runinfo_G0lim0p15.mat')



%% nond a VS a/h (colored by slp) 
figure()
scatter(t_scale_all,a_tot_all,45,'filled')
% hold on 
% scatter(nond2(ind_slp3),a_tot_all(ind_slp3)./(h_all(ind_slp3)),45,'filled')
% hold on 
% scatter(nond2(ind_slp4),a_tot_all(ind_slp4)./(h_all(ind_slp4)),45,'filled')
% hold on 
% plot(Hs_nond_model,a_nond_model,'LineWidth',3.5,'Color','k')
% hold off
xlabel('$\sqrt{h/g} (s)$','Interpreter','latex')
ylabel('a (s)')
%title('$cxt_{fit} = \exp{(-\frac{x}{a})} \cos{(\frac{x}{b} +c)}/\cos{c}$','interpreter','latex')
%title({['a/h = ',num2str(f_a_nond.A),'*Hs/h +',num2str(f_a_nond.B)],[' r^2 = ',num2str(gof_a_nond.rsquare)]},'interpreter','latex')
%legend('slp2','slp3','slp4')
niceplot(18)
xlim([0,Inf])
ylim([0,Inf])

width= 30;
height = 20;
set(gcf,'Units','inches','Position',[0,0,width,height])
set(gcf,'visible','off') 
shg
%fig = gcf;
%fig.Position(4) = fig.Position(4)+100;
%fig.ResizeFcn = @(src,event) set(gca,'Position',[0.1, 0.1, 0.8, 0.8])

%exportgraphics(gcf,'/data1/bliu/figures/paper_figures/nond_abc_vs_nondnum.pdf','ContentType','vector')
%hAx = gca;
%hAx.Position = hAx.Position.*[1 1 1 0.95];


%%
figure()
scatter(nond_Hsh,a_nond,50,'filled')
hold on 
plot(nond_Hsh_model,a_nond_model,'LineWidth',5)
xlabel('$\frac{H_{s}}{h}$','Interpreter','latex')
ylabel('$\frac{a}{\sqrt{h/g}}$','Interpreter','latex')
xlim([0 inf])
ylim([0 inf])
title(['r^2=',num2str(gof_a_nond.rsquare)])
hold off
niceplot(20)

%% nond t vs nond G0
figure()
% subplot(121)
% scatter(x_nond_all,a_nond,50,'filled')
% xlabel('$\frac{x}{x_b}$','Interpreter','latex')
% ylabel('$\frac{a}{\sqrt{h/g}}$','Interpreter','latex')
% niceplot(20)
% 
% subplot(122)
scatter(G0_nond.slp2,a_tot.slp2./t_scale.slp2,50,'filled')
hold on 
scatter(G0_nond.slp3,a_tot.slp3./t_scale.slp3,50,'filled')
hold on 
scatter(G0_nond.slp4,a_tot.slp4./t_scale.slp4,50,'filled')
hold on 
plot(nond_G0_model,a_nond_model_G0,'LineWidth',6)
title(['r^2=',num2str(gof_a_nond_G0.rsquare)])
legend('slp2','slp3','slp4')
xlabel('$\frac{G_0}{G_{0max}}$','Interpreter','latex')
ylabel('$\frac{a}{\sqrt{h/g}}$','Interpreter','latex')
niceplot(20)
hold off 

%%
figure()
scatter(x_nond_all,a_tot_all,50,'filled')
xlabel('$\frac{x}{x_b}$','Interpreter','latex')
ylabel('a')
niceplot(20)