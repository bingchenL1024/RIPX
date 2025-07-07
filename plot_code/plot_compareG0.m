% Bingchen Liu May 20, 2024
% This plot (using oneone) compares magnitude of different choice of G0 

clear
load('/data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat')
%load qa_qc_RIPX_NK.mat
load('/data1/bliu/data/Sky_WBfit_qced') 

%%  
figure()
subplot(2,2,1)
plot_vfrc_scatter_col_slp(A.iS2s,A.wfit2.intw,1)
plot_vfrc_scatter_col_slp(A.iS3s,A.wfit3.intw,2)
plot_vfrc_scatter_col_slp(A.iS4s,A.wfit4.intw,3)
plot_oneone 
ylabel('${\int^{k_{max}} \, S_{WB} \, dk_{y}}$','interpreter','latex')
xlabel('${\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
title('With Cutoff','interpreter','latex')
niceplot_nobold_nomintick(20)

subplot(2,2,2)
plot_vfrc_scatter_col_slp(A.iS2tot,A.wfit2.So,1)
plot_vfrc_scatter_col_slp(A.iS3tot,A.wfit3.So,2)
plot_vfrc_scatter_col_slp(A.iS4tot,A.wfit4.So,3)
plot_oneone 
xlabel('${\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
ylabel('${\int \, S_{WB} \, dk_{y}}$','interpreter','latex')
title('Without Cutoff','interpreter','latex')
niceplot_nobold_nomintick(20)


subplot(2,2,3)
plot_vfrc_scatter_col_slp(A.iS2s,A.iS2tot,1)
plot_vfrc_scatter_col_slp(A.iS3s,A.iS3tot,2)
plot_vfrc_scatter_col_slp(A.iS4s,A.iS4tot,3)
plot_oneone 
xlabel('${\int^{k_{max}} \, S_{ky} \, dk_{y}}$','interpreter','latex')
ylabel('${\int \, S_{ky} \, dk_{y}}$','interpreter','latex')
title('Data Spectra Integral','interpreter','latex')
niceplot_nobold_nomintick(20)


subplot(2,2,4)
plot_vfrc_scatter_col_slp(A.wfit2.intw,A.wfit2.So,1)
plot_vfrc_scatter_col_slp(A.wfit3.intw,A.wfit3.So,2)
plot_vfrc_scatter_col_slp(A.wfit4.intw,A.wfit4.So,3)
plot_oneone 
xlabel('${\int^{k_{max}} \, S_{WB} \, dk_{y}}$','interpreter','latex')
ylabel('${\int \, S_{WB} \, dk_{y}}$','interpreter','latex')
title('Weibull Fit Spectra Integral','interpreter','latex')
niceplot_nobold_nomintick(20)

