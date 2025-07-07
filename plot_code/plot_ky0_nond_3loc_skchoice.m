% Bingchen Liu. May 20, 2024
% Plot ky0 and its scaling


% "K": include directional spread and peak wave number for incoming wave
% spectra
% "A": spectra info from model data and Weibull fit
% "SS": basic model run info from S.mat




clear
%load edited_Wall.mat   
load('/data1/bliu/data/plotready_ky0_nond_3loc.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%
g = 9.81;
%beta = 1.2; % Weibull fit exponent
fs = 12;
sz_loc = [-0.75, -0.5, -0.33];
col2 = [0.3 0.3 1;0.9 0.6 0;1 .2 .2];

%% multiply \beta and Ir
figure()
subplot(2,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h,3)
ylabel('$$k_{y0}\, h$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
niceplot_nobold_nomintick(24)

subplot(2,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h_beta_div,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h_beta_div,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h_beta_div,3)
ylabel('$$k_{y0} \, h \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(24)

subplot(2,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h_Ir_div,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h_Ir_div,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h_Ir_div,3)
ylabel('$$k_{y0}\, h \, \mathrm{Ir}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)

subplot(2,3,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw,3)
ylabel('$$k_{y0} /k_w$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)


subplot(2,3,5)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw_beta_div,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw_beta_div,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw_beta_div,3)
ylabel('$$(k_{y0} /k_w) \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)

subplot(2,3,6)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw_Ir_div,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw_Ir_div,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw_Ir_div,3)
ylabel('$$(k_{y0} /k_w) \, \mathrm{Ir}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)

%% divide by \beta and Ir
figure()
subplot(2,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h,3)
ylabel('$$k_{y0}\, h$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
niceplot_nobold_nomintick(24)

subplot(2,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h_beta_X,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h_beta_X,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h_beta_X,3)
ylabel('$$k_{y0} \, h \, \frac{1}{\beta}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(24)

subplot(2,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_h_Ir_X,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_h_Ir_X,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_h_Ir_X,3)
ylabel('$$k_{y0}\, h \, \frac{1}{\mathrm{Ir}}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)

subplot(2,3,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw,3)
ylabel('$$k_{y0} /k_w$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)


subplot(2,3,5)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw_beta_X,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw_beta_X,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw_beta_X,3)
ylabel('$$(k_{y0} /k_w) \, \frac{1}{\beta}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)

subplot(2,3,6)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./sk_ky0.Sk2_kw_Ir_X,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./sk_ky0.Sk3_kw_Ir_X,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./sk_ky0.Sk4_kw_Ir_X,3)
ylabel('$$(k_{y0} /k_w) \, \frac{1}{\mathrm{Ir}}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.003])
niceplot_nobold_nomintick(24)
