
% Bingchen Liu Nov 22, 2024
% This code is used to plot along-shore wavenumber spectra 
% Correspond to Fig 5 without the Weibull fit,
% NOTE: this code uses specific run, for general spec plot, check 

% data from Master_analysis code --Analysis code needed for calculation of nond G0 -- 10 locaitons 


clear
close all 
% load('/data1/bliu/data/Sky_nond')
% load('/data1/bliu/data/Sky_WBfit_qced') 
% load('/data1/bliu/data/plotready_ky0_nond_3loc.mat')

load('/data1/bliu/data/Sky_binmean_10loc_1p4.mat') %'get_nond_binmean_spectra_10loc'
load('/data1/bliu/data/plotready_ky0_nond_10loc_2025.mat')

%%
xsize=13;ysize=15;
x0=0.145; 
y0=0.095;  
dy=0.1;
xw=0.8; 
yw=0.39;
y1= y0+yw+dy; 
pos = [x0 y0 xw yw; x0 y1 xw yw];


tic

%col = [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880];
fig_fontsize = 10;
legend_size = 8;
subfiglabel_fontsz = 14;
subfiglabel = {'(a)','(b)'};
sub1_linewidth = 0.5;
sub2_linewidth = 0.5;
wb_linewidth  = 1.5;
errbar_col = [1,0,1];

clear col
x_num = length(ky2_nond(1,:,1));
col = cmocean('thermal',x_num);

x_nond = linspace(-0.75,-0.25,5);


figure(1)
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:))

for xloc =x_num:-1:1
    for N =1:24
        h3=loglog(ky_full,squeeze(Sky2_full_runmean(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth);
        hold on 
        loglog(ky_full,squeeze(Sky3_full_runmean(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth)
        hold on 
        loglog(ky_full,squeeze(Sky4_full_runmean(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth)
        hold on 
    end
end 
xline(0.2,'LineWidth',2.5,'Color',[.7 .7 .7],'LineStyle','--')

% clear col
% colormap(cmocean('thermal',x_num))
% cbar = colorbar;
% cbar.Ticks =linspace(-0.75,-0.25,5);
% cbar.Label.FontSize = 16;
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$x / x_{\mathrm{b}}$';

%set(gca,'XTickLabelRotation',0);
hold off 
ax= gca;
%ax.XLabel.Position(2) = -0.069;
ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}}\,(\mathrm{m} \,s^{-4})$','Interpreter','latex')
xlabel('$k_{y} \, (\mathrm{m}^{-1})$','Interpreter','latex');
set(gca,'XTick',[0.001,0.01,0.1])
set(gca,'YTick',[0.00001,0.0001,0.001,0.01,0.1])
set(ax.XLabel,'Position',[0.0255, 1.6617e-06, -3.5])
ylim([0.000005,0.1])
xlim([ky_full(1), 0.4])
grid on 
niceplot_nobold_nomintick(fig_fontsize)     
text(0.00097,0.04,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')



% ===================================================================
%%
col = cmocean('thermal',x_num);

subplot("Position",pos(1,:))
for xloc =length(ky2_nond(1,:,1)):-1:1
    for N =1:24
        h3=loglog(squeeze(ky2_nond(N,xloc,:))',squeeze(Sky2_nond_kym_S0(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth);
        hold on 
        loglog(squeeze(ky3_nond(N,xloc,:))',squeeze(Sky3_nond_kym_S0(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth)
        hold on 
        loglog(squeeze(ky4_nond(N,xloc,:))',squeeze(Sky4_nond_kym_S0(N,xloc,:))','Color',col(xloc,:),'LineWidth',sub2_linewidth)
        hold on 
    end
end 

errorbar(ky_bincenter,Sky_binmean,error_binmean,'s','LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor',errbar_col,'color',errbar_col)
set(gca,'XScale','log','YScale','log')
hold on 
plot(ky_nond,Sky_nond_analyt,'LineWidth',3,'Color','r','LineStyle','--')
hold on 
%scatter(ky_bincenter,Sky_logbinmean,50,'yellow','diamond','filled')
ylabel('$\tilde{S}_{\nabla \times \mathbf{F}_{\mathrm{br}}}$','Interpreter','latex')
%ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, \hat{k}_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
xlabel('$\tilde{k}_{y}$','Interpreter','latex')
%xlim([-inf 10])
%title(['Expo = ',expo_name])
set(gca,'XTick',[0.001,0.01,0.1,1,10])
grid on 
axis([0.03 10 0.01 1])
%cbar.Label.Position = [0,-2.5,0];
clear col
colormap(cmocean('thermal',x_num))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.Label.FontSize = 14;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Position = [0.3, 0.12,0.5,0.02];

%cbar.AxisLocation='out';

%title({['r^2=',num2str(R_square)],['RMSE=',num2str(rmse_bin,2),', RMSE_{logbin}=',num2str(rmse_logbin,2)],['linRMSE_{linmean}=',num2str(rmse_lin,2),', linRMSE_{logbin} =', num2str(rmse_lin_logbin,2)]})
niceplot_nobold_nomintick(fig_fontsize)
text(0.035,0.6,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

hold off


%%

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf -painters '/data1/bliu/vortpaper/fig_sky.pdf'
toc 
close all




%% AGU 2024 version
% for N= 1:24 %original
% loglog(squeeze(ky2_cutoff(N,:,:))',squeeze(Sky2_cutoff_runmean(N,:,:))','LineWidth',sub1_linewidth,'Color',[.5 .5 .5])
% hold on 
% loglog(squeeze(ky3_cutoff(N,:,:))',squeeze(Sky3_cutoff_runmean(N,:,:))','LineWidth',sub1_linewidth,'Color',[.5 .5 .5])
% hold on 
% loglog(squeeze(ky4_cutoff(N,:,:))',squeeze(Sky4_cutoff_runmean(N,:,:))','LineWidth',sub1_linewidth,'Color',[.5 .5 .5])
% end 
% for N= 1:24
% p3=loglog(squeeze(ky2_cutoff(N,3,:))',squeeze(Sky2_cutoff_runmean(N,3,:))','LineWidth',sub1_linewidth,'Color',col(3,:));
% hold on 
% loglog(squeeze(ky3_cutoff(N,3,:))',squeeze(Sky3_cutoff_runmean(N,3,:))','LineWidth',sub1_linewidth,'Color',col(3,:))
% hold on 
% loglog(squeeze(ky4_cutoff(N,3,:))',squeeze(Sky4_cutoff_runmean(N,3,:))','LineWidth',sub1_linewidth,'Color',col(3,:))
% hold on 
% end 
% for N= 1:24
% p1=loglog(squeeze(ky2_cutoff(N,1,:))',squeeze(Sky2_cutoff_runmean(N,1,:))','LineWidth',sub1_linewidth,'Color',col(1,:));
% hold on 
% loglog(squeeze(ky3_cutoff(N,1,:))',squeeze(Sky3_cutoff_runmean(N,1,:))','LineWidth',sub1_linewidth,'Color',col(1,:))
% hold on 
% loglog(squeeze(ky4_cutoff(N,1,:))',squeeze(Sky4_cutoff_runmean(N,1,:))','LineWidth',sub1_linewidth,'Color',col(1,:))
% hold on 
% end 
% for N= 1:24
% p2=loglog(squeeze(ky2_cutoff(N,2,:))',squeeze(Sky2_cutoff_runmean(N,2,:))','LineWidth',sub1_linewidth,'Color',col(2,:));
% hold on 
% loglog(squeeze(ky3_cutoff(N,2,:))',squeeze(Sky3_cutoff_runmean(N,2,:))','LineWidth',sub1_linewidth,'Color',col(2,:))
% hold on 
% loglog(squeeze(ky4_cutoff(N,2,:))',squeeze(Sky4_cutoff_runmean(N,2,:))','LineWidth',sub1_linewidth,'Color',col(2,:))
% hold on 
% end 
% 
% 
% hold off 
% ax= gca;
% %ax.XLabel.Position(2) = -0.069;
% ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}}(s^{-4}/\mathrm{cpm})$','Interpreter','latex')
% xlabel('$k_{y} (\mathrm{cpm})$','Interpreter','latex');
% set(gca,'XTick',[0.01,0.1,1,10])
% set(gca,'YTick',[0.00001,0.0001,0.001,0.01,0.1])
% set(ax.XLabel,'Position',[0.0195, 2.6617e-06, -1])
% ylim([0.00001,0.1])
% xlim([0.0017 0.2])
% grid on 
% niceplot_nobold_nomintick(fig_fontsize)     
% text(0.00185,0.04,subfiglabel{1},'FontSize',subfiglabel_fontsz)
% 
% 
% 
% subplot("Position",pos(1,:))
% 
% for N =1:24
% h3=loglog(squeeze(ky2_nond(N,3,:))',squeeze(Sky2_nond_kym_S0(N,3,:))','Color',col(3,:),'LineWidth',sub2_linewidth);
% hold on 
% loglog(squeeze(ky3_nond(N,3,:))',squeeze(Sky3_nond_kym_S0(N,3,:))','Color',col(3,:),'LineWidth',sub2_linewidth)
% hold on 
% loglog(squeeze(ky4_nond(N,3,:))',squeeze(Sky4_nond_kym_S0(N,3,:))','Color',col(3,:),'LineWidth',sub2_linewidth)
% hold on 
% end
% 
% 
% for N =1:24
% h1=loglog(squeeze(ky2_nond(N,1,:))',squeeze(Sky2_nond_kym_S0(N,1,:))','Color',col(1,:),'LineWidth',sub2_linewidth);
% hold on 
% loglog(squeeze(ky3_nond(N,1,:))',squeeze(Sky3_nond_kym_S0(N,1,:))','Color',col(1,:),'LineWidth',sub2_linewidth)
% hold on 
% loglog(squeeze(ky4_nond(N,1,:))',squeeze(Sky4_nond_kym_S0(N,1,:))','Color',col(1,:),'LineWidth',sub2_linewidth)
% hold on 
% end
% 
% for N =1:24
% h2=loglog(squeeze(ky2_nond(N,2,:))',squeeze(Sky2_nond_kym_S0(N,2,:))','Color',col(2,:),'LineWidth',sub2_linewidth);
% hold on
% loglog(squeeze(ky3_nond(N,2,:))',squeeze(Sky3_nond_kym_S0(N,2,:))','Color',col(2,:),'LineWidth',sub2_linewidth)
% hold on 
% loglog(squeeze(ky4_nond(N,2,:))',squeeze(Sky4_nond_kym_S0(N,2,:))','Color',col(2,:),'LineWidth',sub2_linewidth)
% hold on 
% end
% 
% for N =1:24
% loglog(squeeze(ky2_nond_wb(N,:,:))',squeeze(Sky2_nond_kym_wb_S0(N,:,:))','LineWidth',wb_linewidth,'Color','k','LineStyle','--')
% hold on
% loglog(squeeze(ky3_nond_wb(N,:,:))',squeeze(Sky3_nond_kym_wb_S0(N,:,:))','LineWidth',wb_linewidth,'Color','k','LineStyle','--')
% hold on 
% loglog(squeeze(ky4_nond_wb(N,:,:))',squeeze(Sky4_nond_kym_wb_S0(N,:,:))','LineWidth',wb_linewidth,'Color','k','LineStyle','--')
% hold on 
% end
% 
% 
% ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
% xlabel('$k_{y}/k_{y0}$','Interpreter','latex')
% %xlim([-inf 10])
% set(gca,'XTick',[0.001,0.01,0.1,1,10])
% grid on 
% axis([0.03 10 0.01 1])
% niceplot_nobold_nomintick(fig_fontsize)
% %leg=legend([h1,h2,h3],{'$-0.75 \, L_{sz}$','$-0.5 \, L_{sz}$','$-0.33 \, L_{sz}$'},'FontSize',legend_size,'Location','northeast','Interpreter','latex');
% %leg.ItemTokenSize(1) = 12;
% %leg.ItemTokenSize(2) = 2;
% leg=legend([h1,h2,h3],{'$-0.75 \, L_{sz}$','$-0.5 \, L_{sz}$','$-0.33 \, L_{sz}$'},'FontSize',legend_size,'Interpreter','latex');
% leg.ItemTokenSize(1) = 12;
% leg.Position =[ 0.82, 0.38, 0.04, 0.04];  
% text(0.035,0.6,subfiglabel{2},'FontSize',subfiglabel_fontsz)
% hold off
% 
% toc
% set(gcf,'PaperUnit','centimeters')
% set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])
% 
% print -dpdf -painters '/data1/bliu/vortpaper/fig_sky.pdf'
% toc 
% close all


