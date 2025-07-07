% Bingchen Liu Nov 29,2024
% This code plot the nond G0 and nond ky0 scaling vs directional spread 
% modified from 'plot_scal_nond' --> without saving actual figure 



clearvars -except expo


G0=load('/data1/bliu/data/plotready_G0_nond_3loc_2025.mat');
ky0=load('/data1/bliu/data/plotready_ky0_nond_3loc_2025.mat');



% xsize=8;ysize=11;
% x0=0.23;
% y0=0.117;
% dy=0.05;
% xw=0.7;
% yw=0.4108;
% y1= y0+yw+dy;
% pos = [x0 y0 xw yw; x0 y1 xw yw];



axfont_sz = 25;
fig_fontsz= 25;
subfiglabel_fontsz = 15;
subfiglabel = {'(a)','(b)'};


figure()

subplot(211);
plot(G0.ds_model,G0.G0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
plot_vfrc_scatter_col_slp(repmat(G0.ds.sigtb2,[3 1])',G0.G0_nond.slp2,1)
plot_vfrc_scatter_col_slp(repmat(G0.ds.sigtb3,[3 1])',G0.G0_nond.slp3,2)
plot_vfrc_scatter_col_slp(repmat(G0.ds.sigtb4,[3 1])',G0.G0_nond.slp4,3)
%xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
%title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$\frac{G_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}_{\infty}$','interpreter','latex','fontsize',axfont_sz);
ax= gca;
set(gca,'XTickLabel',[])
%ax.YLim = [0 inf];
title(['expo =',num2str(expo),', r^2 =',num2str(G0.gof.rsquare)])
xlim([0,15])
ylim([0, 14])
hold off

subplot(212);
plot(ky0.ds_model,ky0.ky0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
plot_vfrc_scatter_col_slp(repmat(ky0.ds.sigtb2,[3 1])',ky0.ky0_nond.slp2,1)
plot_vfrc_scatter_col_slp(repmat(ky0.ds.sigtb3,[3 1])',ky0.ky0_nond.slp3,2)
plot_vfrc_scatter_col_slp(repmat(ky0.ds.sigtb4,[3 1])',ky0.ky0_nond.slp4,3)
%title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(fig_fontsz)
title(['r^2 =',num2str(ky0.gof.rsquare)])
ylabel('$$k_{y0} \, h \, \beta$$','interpreter','latex','fontsize',axfont_sz);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex','fontsize',axfont_sz);
hold off


