% Bingchen Liu. May 20, 2024
% 8 panels plot 
% This code plot the first version of nondF with 3 crossshore location and
% it also includes plot of dimensional G0 with different definition of G0
% for visual 


clear
close all
%addpath(genpath('/home/ssuanda/matlab/')) %modified by BL 
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
addpath(genpath('/data1/bliu'))
addpath(genpath('/home/ffeddersen/matlab'))


load('/data1/bliu/data/qa_qc_RIPX_NK.mat') %load qa_qc_RIPX_NK.mat
%load edited_Wall.mat   
%SS = load('/data1/bliu/data/S0scale_10');
load('/data1/bliu/data/plotready_G0_nond_3loc')

%% calculating scaling 
% %for plotting
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% g = 9.81;
% beta = 1.2; % Weibull fit exponent
% fs = 12;
% col2 = [0.3 0.3 1;0.9 0.6 0;1 .2 .2];
% sz_loc = [-0.75, -0.5, -0.33];
% 
% %%%%%%%%%%%%%%% scaling 
% Sk2 = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5));
% Sk3 = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5));
% Sk4 = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5));
% 
% Sk22 = SS.Dw2./(SS.h2.*((g.*(SS.h2.^3)).^.5));
% Sk32 = SS.Dw3./(SS.h3.*((g.*(SS.h3.^3)).^.5));
% Sk42 = SS.Dw4./(SS.h4.*((g.*(SS.h4.^3)).^.5));
%  
% Sk23 = SS.Dw2./((SS.h2.*((g.*(SS.h2.^3)).^.5)).*Ir2); %include Iribarren number
% Sk33 = SS.Dw3./((SS.h3.*((g.*(SS.h3.^3)).^.5)).*Ir3); %include Iribarren number
% Sk43 = SS.Dw4./((SS.h4.*((g.*(SS.h4.^3)).^.5)).*Ir4); %include Iribarren number
% 
% % Sk23 = (A.kym2.*SS.Dw2)./(((g.*(SS.h2.^3)).^.5));
% % Sk33 = (A.kym3.*SS.Dw3)./(((g.*(SS.h3.^3)).^.5));
% % Sk43 = (A.kym4.*SS.Dw4)./(((g.*(SS.h4.^3)).^.5));
% 

%% plot: dimensional curl
% figure()
% plot_vfrc_scatter_col_slp(repmat(K.sigtb2,[3 1])',sqrt(A.wfit2.intw),1)
% plot_vfrc_scatter_col_slp(repmat(K.sigtb3,[3 1])',sqrt(A.wfit3.intw),2)
% plot_vfrc_scatter_col_slp(repmat(K.sigtb4,[3 1])',sqrt(A.wfit4.intw),3)
% ylabel('$\sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex','fontsize',16);
% xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
% title('Dimensional $\nabla \times \mathbf{F}_{\mathrm{br}}$','interpreter','latex')
% niceplot_nobold_nomintick(24)
% xlim([0,15])


%% 8 panels of nondF w/ and w/o \beta, 1/h as inverse length scale

figure()
subplot(3,4,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.iS2tot)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.iS3tot)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.iS4tot)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.iS2s)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.iS3s)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.iS4s)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.So)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.So)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.So)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_h,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_h,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_h,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

%%%%%%%%%%%%%%%%%% Multiply by \beta

subplot(3,4,5)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2tot)./Sk.Sk2_h).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3tot)./Sk.Sk3_h).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4tot)./Sk.Sk4_h).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h \beta}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,6)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2s)./Sk.Sk2_h).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3s)./Sk.Sk3_h).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4s)./Sk.Sk4_h).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h \beta}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,7)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.So)./Sk.Sk2_h).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.So)./Sk.Sk3_h).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.So)./Sk.Sk4_h).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h \beta}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int  S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(3,4,8)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h \beta}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


%%%%%%%%%%%%%%%%%% Multiply by Iribarren number

subplot(3,4,9)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2tot)./Sk.Sk2_h_Ir).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3tot)./Sk.Sk3_h_Ir).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4tot)./Sk.Sk4_h_Ir).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w} \, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,10)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2s)./Sk.Sk2_h_Ir).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3s)./Sk.Sk3_h_Ir).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4s)./Sk.Sk4_h_Ir).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}\, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,11)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.So)./Sk.Sk2_h_Ir).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.So)./Sk.Sk3_h_Ir).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.So)./Sk.Sk4_h_Ir).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w} \, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int  S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(3,4,12)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_h_Ir),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_h_Ir),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_h_Ir),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h }{D_\mathrm{w}\, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])





%% 8 panels of nondF w/ and w/o \beta, kw as inverse length scale

figure()
subplot(3,4,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.iS2tot)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.iS3tot)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.iS4tot)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.iS2s)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.iS3s)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.iS4s)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,3)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.So)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.So)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.So)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,4)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',sqrt(A.wfit2.intw)./Sk.Sk2_k,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',sqrt(A.wfit3.intw)./Sk.Sk3_k,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',sqrt(A.wfit4.intw)./Sk.Sk4_k,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

%%%%%%%%%%%%%%%%%% Multiply by \beta

subplot(3,4,5)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2tot)./Sk.Sk2_k).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3tot)./Sk.Sk3_k).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4tot)./Sk.Sk4_k).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} \beta}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,6)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2s)./Sk.Sk2_k).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3s)./Sk.Sk3_k).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4s)./Sk.Sk4_k).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} \beta}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,7)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.So)./Sk.Sk2_k).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.So)./Sk.Sk3_k).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.So)./Sk.Sk4_k).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} \beta}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int  S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])
%ylim([0,10])


subplot(3,4,8)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k).*0.02,1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k).*0.03,2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k).*0.04,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} \beta}{D_\mathrm{w} \, k_{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])
%ylim([0,10])

%%%%%%%%%%%%%%%%%% Multiply by Ir

subplot(3,4,9)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2tot)./Sk.Sk2_k_Ir),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3tot)./Sk.Sk3_k_Ir),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4tot)./Sk.Sk4_k_Ir),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}\, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,10)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.iS2s)./Sk.Sk2_k_Ir),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.iS3s)./Sk.Sk3_k_Ir),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.iS4s)./Sk.Sk4_k_Ir),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}\, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])

subplot(3,4,11)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.So)./Sk.Sk2_k_Ir),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.So)./Sk.Sk3_k_Ir),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.So)./Sk.Sk4_k_Ir),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w}\, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int  S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


subplot(3,4,12)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb2,[3 1])',(sqrt(A.wfit2.intw)./Sk.Sk2_k_Ir),1)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb3,[3 1])',(sqrt(A.wfit3.intw)./Sk.Sk3_k_Ir),2)
plot_vfrc_scatter_col_slp(repmat(ds.sigtb4,[3 1])',(sqrt(A.wfit4.intw)./Sk.Sk4_k_Ir),3)
ylabel('$\frac{G_{0} (gh^3)^{1/2}}{D_\mathrm{w} \, k_{w} \, \mathrm{Ir}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])




