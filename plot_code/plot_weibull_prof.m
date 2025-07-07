% Bingchen Liu, March/July 2024
% modified from Ata's original code 'plot_all_run_weibul'
% This code plots and test different scaling as well as Weibull fit
% It will plot above variables as a function of cross-shore location and
% see how good the scaling is at diff location
% It will loop for all the good runs and create/save figue for each run

clf;close all;clear;
warning off 
load('/data1/bliu/data/Profs_fit_to0p2_BL');

for i = 2:4
slp = i; % Change this to a different slope

for N = 1:24;
eval(sprintf('A = Al_s%1.0f.a%1.0f;',slp,N))
fname = sprintf('sl%1.0f_hs%1.1fm_tp%1.0fs_ds%2.0fd',slp,A.Hs(1),A.inp.tp,A.inp.dsp);
stit = sprintf('sl=%1.0f hs=%1.1fm tp=%1.0fs ds=%2.0fd',slp,A.Hs(1),A.inp.tp,A.inp.dsp);
Sk1 = (A.Dw)./(((9.81.*(A.h)).^.5.*A.hb.^2.*A.beta)); % scaling with hb (for hb^2)
%Sk3 = (A.kp.^2.*A.Dw)./(((9.81.*(A.h)).^.5));
kh = A.kp.*A.h;

xl = min(A.X);fs = 22;
set(gcf,'color','w','position',[7 36 515 664])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dim profile plot
figure(1)
subplot(2,2,1)
plot(A.X,-A.h,'k','linewidth',2)
hold on
plot(A.X3locs,-A.h3locs,'or','markerfacecolor','r')
xlim([xl 0]);grid on;ylabel('h (m)')
title(stit)
niceplot_nobold_nomintick(fs);

subplot(2,2,3)
plot(A.X,A.raw.kyo,'k','linewidth',2)
hold on
plot(A.X,A.web.kyo,'r','linewidth',2)
plot(A.X,A.kp,'b','linewidth',2)
legend('Model','WB Fit','kp','location','northwest')
xlim([xl 0]);grid on;ylabel('ky0 (m^{-1})')
niceplot_nobold_nomintick(fs);

subplot(2,2,2)
plot(A.X,A.web.lsk,'k','linewidth',2)
hold on
plot(A.X,A.web.wsk,'r','linewidth',2)
legend('vsk','wsk','location','northwest')
xlim([xl 0]);ylim([0 1]);grid on;ylabel('Skill)')
niceplot_nobold_nomintick(fs);

subplot(2,2,4)
yyaxis left;
plot(A.X,sqrt(A.raw.is),'k','linewidth',2)
hold on
plot(A.X,sqrt(A.web.is),'r','linewidth',2)
ax= gca;
ax.YColor = 'k';
ylabel('$\sqrt{\int^{full} \, S_{ky} \, dk_{y}} \, (s^{-2})$','interpreter','latex');

yyaxis right
plot(A.X, Sk1,'b','linewidth',2)
ylabel('$\frac{D_w}{c_{l} h_{b}^2} \frac{1}{Ir} \, (s^{-2})$','interpreter','latex');
ax= gca;
ax.YColor = 'b';
legend('Model','WB Fit','Scaling','location','northwest')
xlim([xl 0]);grid on;
niceplot_nobold_nomintick(fs);

width= 30;
height = 15;
set(gcf,'Units','inches','Position',[0,0,width,height])
%set(gcf,'visible','off')
saveas(gcf,['/data1/bliu/figures/Weibull_profs/dimensional/',fname,'.png'])
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 nondim profile plot
figure(2)
subplot(3,1,1)
plot(A.X,-A.h,'k','linewidth',2)
hold on
plot(A.X3locs,-A.h3locs,'or','markerfacecolor','r')
xlim([xl 0]);grid on;ylabel('h (m)')
title(stit)
niceplot_nobold_nomintick(fs);

subplot(3,1,2)
plot(A.X,(A.raw.kyo).*A.h.*A.beta,'k','linewidth',2)
hold on
plot(A.X,(A.web.kyo).*A.h.*A.beta,'r','linewidth',2)
legend('Model','WB Fit','location','northwest')
xlim([xl 0]);grid on;ylabel('$k_{y0} h$','interpreter','latex')
niceplot_nobold_nomintick(fs);

subplot(3,1,3)
plot(A.X,(sqrt(A.raw.is)./(Sk1)).*A.Ir,'k','linewidth',2)
hold on
plot(A.X,(sqrt(A.web.is)./(Sk1)).*A.Ir,'r','linewidth',2)
legend('Model','WB Fit','location','northwest')
xlim([xl 0]);grid on;
ylabel('$\frac{\sqrt{\int^{full} \, S_{ky} \, dk_{y}}}{D_{w}/(c_{l} h_{b}^2)}$','interpreter','latex')
niceplot_nobold_nomintick(fs);

width= 30;
height = 15;
set(gcf,'Units','inches','Position',[0,0,width,height])
set(gcf,'visible','off')
saveas(gcf,['/data1/bliu/figures/Weibull_profs/nondimensional/',fname,'.png'])

clf 
close all 
end 
end