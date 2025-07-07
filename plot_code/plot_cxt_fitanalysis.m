% Bingchen Liu, Nov 4, 2024
% This code analyze the a, b, c free parameters and try to develop their
% scaling
% only include the good ones for paper 
clear
%load('/data1/bliu/data/cxt_x_fitanalysis')
load('/data1/bliu/data/cxt_nointerp_fitanalysis')
%load('/data1/bliu/data/cxt_x_maxvar_fitanalysis')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_runinfo')



%% nond a, b, and c VS a/h (colored by slp) -------with linear fit 
figure()
subplot(1,3,1)
scatter(nond2(ind_slp2),a_tot_all(ind_slp2)./(h_all(ind_slp2)),45,'filled')
hold on 
scatter(nond2(ind_slp3),a_tot_all(ind_slp3)./(h_all(ind_slp3)),45,'filled')
hold on 
scatter(nond2(ind_slp4),a_tot_all(ind_slp4)./(h_all(ind_slp4)),45,'filled')
hold on 
plot(Hs_nond_model,a_nond_model,'LineWidth',3.5,'Color','k')
hold off
xlabel('Hs/h')
ylabel('a/h (decay rate)')
%title('$cxt_{fit} = \exp{(-\frac{x}{a})} \cos{(\frac{x}{b} +c)}/\cos{c}$','interpreter','latex')
title({['a/h = ',num2str(f_a_nond.A),'*Hs/h +',num2str(f_a_nond.B)],[' r^2 = ',num2str(gof_a_nond.rsquare)]},'interpreter','latex')
legend('slp2','slp3','slp4')
niceplot(18)
xlim([0,Inf])
ylim([0,Inf])

subplot(1,3,2)
scatter(nond2(ind_slp2),b_tot_all(ind_slp2)./(h_all(ind_slp2)),45,'filled')
hold on 
scatter(nond2(ind_slp3),b_tot_all(ind_slp3)./(h_all(ind_slp3)),45,'filled')
hold on 
scatter(nond2(ind_slp4),b_tot_all(ind_slp4)./(h_all(ind_slp4)),45,'filled')
hold on 
plot(Hs_nond_model,b_nond_model,'LineWidth',3.5,'Color','k')
hold off 
xlabel('Hs/h')
ylabel('b/h (decay rate)')
%title('$cxt_{fit} = \exp{(-\frac{x}{a})} \cos{(\frac{x}{b} +c)}/\cos{c}$','interpreter','latex')
title({['b/h = ',num2str(f_b_nond.A),'*Hs/h +',num2str(f_b_nond.B)],[' r^2 = ',num2str(gof_b_nond.rsquare)]},'interpreter','latex')
legend('slp2','slp3','slp4')
niceplot(18)
xlim([0,Inf])
ylim([0,Inf])

subplot(1,3,3)
scatter(nond2(ind_slp2),c_tot_all(ind_slp2),45,'filled')
hold on 
scatter(nond2(ind_slp3),c_tot_all(ind_slp3),45,'filled')
hold on 
scatter(nond2(ind_slp4),c_tot_all(ind_slp4),45,'filled')
hold on 
plot(Hs_nond_model,c_nond_model,'LineWidth',3.5,'Color','k')
xlabel('Hs/h')
ylabel('c (phase shift)')
legend('slp2','slp3','slp4')
title({['b = ',num2str(f_c_nond.A),'*Hs/h +',num2str(f_c_nond.B)],[' r^2 = ',num2str(gof_c_nond.rsquare)]},'interpreter','latex')
niceplot(18)
xlim([0,Inf])
ylim([0,Inf])

width= 30;
height = 5;
set(gcf,'Units','inches','Position',[0,0,width,height])
set(gcf,'visible','off') 

%fig = gcf;
%fig.Position(4) = fig.Position(4)+100;
%fig.ResizeFcn = @(src,event) set(gca,'Position',[0.1, 0.1, 0.8, 0.8])

%exportgraphics(gcf,'/data1/bliu/figures/paper_figures/nond_abc_vs_nondnum.pdf','ContentType','vector')
%hAx = gca;
%hAx.Position = hAx.Position.*[1 1 1 0.95];
%% a and b with linear fit 

figure()
subplot(1,2,1)
scatter(h_all,a_tot_all,45,'filled')
hold on 
plot(h_model,a_model,'LineWidth',5)
hold off
legend('Data','Linear Fit')
title({['a = ',num2str(f_a.A),'*h+',num2str(f_a.B)],[' r^2 = ',num2str(gof_a.rsquare)]},'interpreter','latex')
xlabel('h(m)')
ylabel('Fit Parameter a (m)')
ylim([0,Inf])
niceplot(18)

subplot(1,2,2)
scatter(h_all,b_tot_all,45,'filled')
hold on 
plot(h_model,b_model,'LineWidth',5)
hold off
legend('Data','Linear Fit')
title({['b = ',num2str(f_b.A),'*h+',num2str(f_b.B)],[' r^2 = ',num2str(gof_b.rsquare)]},'interpreter','latex')
xlabel('h(m)')
ylabel('Fit Parameter b (m)')
niceplot(18)



