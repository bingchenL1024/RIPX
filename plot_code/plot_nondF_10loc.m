% Bingchen Liu. May 16, 2024
% This code plot nondimF at different cross-shore location (include more than 3 locations)
% 



clear

%% Prep 
addpath(genpath('/home/ssuanda/matlab/')) %modified by BL 
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
addpath(genpath('/data1/bliu'))
addpath(genpath('/home/ffeddersen/matlab'))


load('/data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat')
%load qa_qc_RIPX_NK.mat
load edited_Wall.mat   
SS = load('/data1/bliu/data/S0scale_10');

%for plotting
fs = 12;
col2 = [0.3 0.3 1;0.9 0.6 0;1 .2 .2];

%for plotting surf zone cross shore location 
sz_loc = [-0.75, -0.5, -0.33];


%collect section in "plosdft_falk_figs.m"
% X = [K.kl2pm(:)' K.kl3pm(:)' K.kl4pm(:)'];
% X2 = [1./K.h2(:)' 1./K.h3(:)' 1./K.h4(:)'];
% 
% Y = [A.kym2(:)' A.kym3(:)' A.kym4(:)'];
% Y2 = [A.wkym2(:)' A.wkym3(:)' A.wkym4(:)'];
% 
% [reg_est, Sfit, B, sse] = linreg(X,Y,0);
% [reg_est, Sfit2, B, sse] = linreg(X2,Y,0);
% [reg_est, Sfit3, B, sse] = linreg(X,Y2,0);
% [reg_est, Sfit4, B, sse] = linreg(X2,Y2,0);
%%%%%%%%%%%%%%%%%%%%%%%%%
g = 9.81;
beta = 1.2; % Weibull fit exponent
% Try 2 scales:
Scale2 = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5));
Scale3 = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5));
Scale4 = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5));

Scale22 = SS.Dw2./(SS.h2.*((g.*(SS.h2.^3)).^.5));
Scale32 = SS.Dw3./(SS.h3.*((g.*(SS.h3.^3)).^.5));
Scale42 = SS.Dw4./(SS.h4.*((g.*(SS.h4.^3)).^.5));


nondF2 = SS.curlFstd2./Scale22;
nondF3 = SS.curlFstd3./Scale32;
nondF4 = SS.curlFstd4./Scale42;


[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

% X = [Scale2(:)' Scale3(:)' Scale4(:)'];
% X2 = [Scale22(:)' Scale32(:)' Scale42(:)'];
% 
% Y = sqrt([A.iS2s(:)' A.iS3s(:)' A.iS4s(:)']);
% Y2 = sqrt([A.wfit2.So(:)' A.wfit3.So(:)' A.wfit4.So(:)']);
% Y3 = sqrt([A.wfit2.intw(:)' A.wfit3.intw(:)' A.wfit4.intw(:)']);
% 
% [reg_est, Sfit, B, sse] = linreg(X,Y,0);
% [reg_est, Sfit2, B, sse] = linreg(X2,Y,0);
% [reg_est, Sfit3, B, sse] = linreg(X,Y3,0);
% [reg_est, Sfit4, B, sse] = linreg(X2,Y3,0);


%% plot dimensional G0 for visual 
figure()
plot_vfrc_scatter_col_slp(repmat(K.sigtb2,[3 1])',sqrt(A.iS2tot),1)
plot_vfrc_scatter_col_slp(repmat(K.sigtb3,[3 1])',sqrt(A.iS3tot),2)
plot_vfrc_scatter_col_slp(repmat(K.sigtb4,[3 1])',sqrt(A.iS4tot),3)
ylabel('Integrated total var of DATA spectra ($S^-2$)','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Dimensional $\nabla \times \mathbf{F}_{\mathrm{br}}$','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])


%% plot nondimF

figure()
plot_vfrc_scatter_noshape(repmat(SS.sigtb2,[1 10]),SS.curlFstd2./Scale22,1)
hold on 
plot_vfrc_scatter_noshape(repmat(SS.sigtb3,[1 10]),SS.curlFstd3./Scale32,2)
hold on 
plot_vfrc_scatter_noshape(repmat(SS.sigtb4,[1 10]),SS.curlFstd4./Scale42,3)
ylabel('$\frac{G_{0} (gh^3)^{1/2} h}{D_\mathrm{w}}$','interpreter','latex','fontsize',16);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
title('Scaling standard deviation of curlfF_std','interpreter','latex')
niceplot_nobold_nomintick(24)
xlim([0,15])
hold off





