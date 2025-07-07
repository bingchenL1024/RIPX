% Bingchen Liu, June 17, 2024
% This code investigate the effect of choice of different scaling on nond G0

% "K": include directional spread and peak wave number for incoming wave
% spectra
% "A": spectra info from model data and Weibull fit
% "SS": basic model run info from S.mat

clear
close all
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
addpath(genpath('/data1/bliu'))
addpath(genpath('/home/ffeddersen/matlab'))


load('/data1/bliu/data/qa_qc_RIPX_NK.mat') %load qa_qc_RIPX_NK.mat
load('/data1/bliu/data/plotready_G0_nond_3loc')

%  Plot  different scaling length scale (h and k_w) and different multiplication factor (beta and Ir and gamma)
%% best one so far 
figure()
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div),3)
ylabel('$\frac{G_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
ax.YLim = [0 inf];
xlim([0,15])

%% compare with the best 
figure()
subplot(1,2,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div),3)
ylabel('$\frac{G_{0} c_{linear} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
%ax.YLim = [0 inf];
ylim([0 inf])
xlim([0,15])

subplot(1,2,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div_solit),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div_solit),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div_solit),3)
ylabel('$\frac{G_{0} c_{solit} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])
%% nond G0 * \beta or Ir

%%%%%%%%%%%%%%%%%%%%% Using inverse length scale hb, h at breaking
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_hb,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_hb,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_hb,3)
ylabel('$\frac{G_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_beta_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_beta_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_beta_div),3)
ylabel('$\frac{G_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \beta$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(1,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div),3)
ylabel('$\frac{G_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
ax.YLim = [0 inf];
xlim([0,15])


%%%%%%%%%%%%%%%%%%%%%%%% using solitary wave velocity
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_c_solit,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_c_solit,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_c_solit,3)
ylabel('$\frac{G_{0}}{D_\mathrm{w}/ (c h h_{b})}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$', '$\\ c = \sqrt{g(h + H_{s})}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_c_solit_beta_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_c_solit_beta_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_c_solit_beta_div),3)
ylabel('$\frac{G_{0}}{D_\mathrm{w}/ (c h h_{b})} \beta$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$', '$\\ c = \sqrt{g(h + H_{s})}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(1,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_c_solit_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_c_solit_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_c_solit_Ir_div),3)
ylabel('$\frac{G_{0}}{D_\mathrm{w}/ (c h h_{b})} Ir$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$', '$\\ c = \sqrt{g(h + H_{s})}$','interpreter','latex')
niceplot_nobold_nomintick(24)
ax.YLim = [0 inf];
xlim([0,15])



%%%%%%%%%%%%%%%%%%%%%%%%% Using inverse length scale 1/h
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_beta_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_beta_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_beta_div),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}} \, \beta$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(1,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_Ir_div),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h }{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])



%%%%%%%%%%%%%%%%%%%%% Using inverse length scale kw
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_beta_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_beta_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_beta_div),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}} \, \beta$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(1,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_Ir_div),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])



%% nond G0 divide \beta or Ir
figure()
subplot(2,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(2,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_beta_X),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_beta_X),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_beta_X),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}} \, \frac{1}{\beta}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(2,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_Ir_X),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_Ir_X),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_Ir_X),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h }{D_\mathrm{w}} \, \frac{1}{\mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])



%%%%%%%%%%%%%%%%%%%%% Using inverse length scale kw
subplot(2,3,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(2,3,5)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_beta_X),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_beta_X),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_beta_X),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}} \, \frac{1}{\beta}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(2,3,6)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_Ir_X),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_Ir_X),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_Ir_X),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}} \, \frac{1}{\mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


%% include kh to see how it scales nond G0
figure()
plot_vfrc_scatter_col_szloc((kh.slp2).^2,(sqrt(A.wfit2.intw)./Sk.Sk2_h),1)
hold on 
plot_vfrc_scatter_col_szloc((kh.slp3).^2,(sqrt(A.wfit3.intw)./Sk.Sk3_h),2)
hold on
plot_vfrc_scatter_col_szloc((kh.slp4).^2,(sqrt(A.wfit4.intw)./Sk.Sk4_h),3)
hold off 
xlabel('$(kh)^2$','Interpreter','latex')
ylabel('$\frac{G_{0} (gh^3)^{1/2} \,h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
niceplot(25)

figure()
plot_vfrc_scatter_col_szloc((kh.slp2).^1,(sqrt(A.wfit2.intw)./Sk.Sk2_k),1)
hold on 
plot_vfrc_scatter_col_szloc((kh.slp3).^1,(sqrt(A.wfit3.intw)./Sk.Sk3_k),2)
hold on
plot_vfrc_scatter_col_szloc((kh.slp4).^1,(sqrt(A.wfit4.intw)./Sk.Sk4_k),3)
hold off 
xlabel('$(kh)$','Interpreter','latex')
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
niceplot(25)

%% nond G0 with kh VS dirsp

figure()
subplot(121)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_kh),1)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_kh),2)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_kh),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w} (kh)^4}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');

subplot(122)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_kh),1)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_kh),2)
plot_vfrc_scatter_col_szloc(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_kh),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w} (kh)^5}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');


%% nond G0 including 

%%%%%%%%%%%%%%%%%%%%% Using inverse length scale kw
figure()
subplot(1,3,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div),3)
ylabel('$\frac{G_{0} c_{linear} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
%ax.YLim = [0 inf];
ylim([0 inf])
xlim([0,15])

subplot(1,3,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div_solit),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div_solit),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div_solit),3)
ylabel('$\frac{G_{0} c_{solit} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(1,3,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div_solit_gamma),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div_solit_gamma),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div_solit_gamma),3)
ylabel('$\frac{G_{0} c_{solit} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir} \gamma^3$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])




