% Bingchen Liu. May 20, 2024
% Plot ky0 using Model data and Weibull fit data 

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

%% Dimensional ky0 with directional spread
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.kym2,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.kym3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.kym4,3)
ylabel('$k_{y0\mathrm{\_data}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$k_{y0}$ from Model Data','interpreter','latex')
ylim([0,0.12])
niceplot_nobold_nomintick(22)

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4,3)
ylabel('$$k_{y0\mathrm{\_Weibull}}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$k_{y0}$ from Weibull Fit','interpreter','latex')
ylim([0,0.12])
niceplot_nobold_nomintick(22)

subplot(1,3,3)
plot_vfrc_scatter_col_slp(A.kym2,A.wkym2,1)
plot_vfrc_scatter_col_slp(A.kym3,A.wkym3,2)
plot_vfrc_scatter_col_slp(A.kym4,A.wkym4,3)
plot_oneone
ylabel('$k_{y0\mathrm{\_Weibull}}$','interpreter','latex','fontsize',16);
xlabel('$k_{y0\mathrm{\_data}}$','interpreter','latex');
title('$k_{y0}$ from Data VS Weibull Fit','interpreter','latex')
niceplot_nobold_nomintick(22)

%% 4 panel of nond ky0 using 1/h
figure()
subplot(2,2,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.kym2./scale2_ky0,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.kym3./scale3_ky0,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.kym4./scale4_ky0,3)
ylabel('$$k_{y0\mathrm{\_dat}}\, h$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Model Data $k_{y0}$','interpreter','latex')
xlim([0,15])
niceplot_nobold_nomintick(24)

subplot(2,2,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./scale2_ky0,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./scale3_ky0,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./scale4_ky0,3)
ylabel('$$k_{y0\mathrm{\_wb}} \, h$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,0.1])
niceplot_nobold_nomintick(24)

subplot(2,2,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(A.kym2./scale2_ky0).*beta(2),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(A.kym3./scale3_ky0).*beta(3),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(A.kym4./scale4_ky0).*beta(4),3)
ylabel('$$k_{y0\mathrm{\_dat}}\, h \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Model Data $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,0.003])
niceplot_nobold_nomintick(24)

subplot(2,2,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(A.wkym2./scale2_ky0).*beta(2),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(A.wkym3./scale3_ky0).*beta(3),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(A.wkym4./scale4_ky0).*beta(4),3)
ylabel('$$k_{y0\mathrm{\_wb}} \,h \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,0.003])
niceplot_nobold_nomintick(24)

%% 4 panel of nond ky0 using kw
figure()
subplot(2,2,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.kym2./scale22_ky0,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.kym3./scale32_ky0,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.kym4./scale42_ky0,3)
ylabel('$$\frac{k_{y0\mathrm{\_dat}}}{k_{w}}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Model Data $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,3])
niceplot_nobold_nomintick(24)

subplot(2,2,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',A.wkym2./scale22_ky0,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',A.wkym3./scale32_ky0,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',A.wkym4./scale42_ky0,3)
ylabel('$$\frac{k_{y0\mathrm{\_wb}}}{k_{w}}$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,3])
niceplot_nobold_nomintick(24)

subplot(2,2,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(A.kym2./scale22_ky0).*beta(2),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(A.kym3./scale32_ky0).*beta(3),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(A.kym4./scale42_ky0).*beta(4),3)
ylabel('$$\frac{k_{y0\mathrm{\_dat}}}{k_{w}} \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Model Data $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,0.1])
niceplot_nobold_nomintick(24)

subplot(2,2,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(A.wkym2./scale22_ky0).*beta(2),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(A.wkym3./scale32_ky0).*beta(3),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(A.wkym4./scale42_ky0).*beta(4),3)
ylabel('$$\frac{k_{y0\mathrm{\_wb}}}{k_{w}} \, \beta$$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
ylim([0,0.1])
niceplot_nobold_nomintick(24)
