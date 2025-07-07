% Bingchen Liu Nov 27, 2024
% This code plot the one to one comparison for G0 and ky0 between model
% data and WB fit

clear
close all
%load('/data1/bliu/data/Sky_nond')
%load('/data1/bliu/data/Sky_WBfit_qced') 

load('data1/bliu/data/Sky_nond_10loc_1p4.mat') %'get_nond_spectra_10loc'
%load('/data1/bliu/data/Sky_binmean_10loc_1p4')


rmse_G0=stats.rmse_G0 % in 'get_weibul_qced_10loc'
rmse_ky0= stats.rmse_ky0
rsquare_G0= stats.rsquare_G0
rsquare_ky0= stats.rsquare_ky0

xsize=16;ysize=8;
x0=0.11; 
y0=0.18;  
dx=0.15;
xw=0.36; yw=xw*(xsize/ysize);
x1= x0+xw+dx; 
pos = [x0 y0 xw yw; x1 y0 xw yw];

fig_fontsz= 10;
subfiglabel_fontsz = 15;
subfiglabel = {'(a)','(b)'};

markersz= 20;



% =============================================

x_num = length(A.iS2s(1,:));
x_nond = linspace(-0.75,-0.25,5);

clear col
%col = cmocean('thermal',x_num);
col_single = [.3 .3 .3];
col = repmat(col_single,5,1);

figure()
h1=subplot("Position",pos(1,:));
for xind = 1:x_num
scatter((A.iS2s(:,xind)).^(0.5),(A.wfit2.intw(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
scatter((A.iS3s(:,xind)).^(0.5),(A.wfit3.intw(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
scatter((A.iS4s(:,xind)).^(0.5),(A.wfit4.intw(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
end 
plot_oneone 
niceplot_nobold(fig_fontsz)
hold off 
ylabel('$\hat{G}_{0} \, (s^{-2})$','Interpreter','latex') %'${\int_{0.001}^{0.2} \, S_{WB} \, dk_{y}} (s^{-4})$','interpreter','latex')
xlabel('$G_{0} \, (s^{-2})$','Interpreter','latex') %'${\int_{0.001}^{0.2} \, S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, dk_{y}} (s^{-4})$','interpreter','latex')
%title(['$r^2=$',num2str(rsquare_G0,3)],'Interpreter','latex')
xlim([0 0.08])
ylim([0 0.08])
ax=gca;
%set(ax.XLabel,'Position',[0.00225,-6.5192e-04,-1])
xtick(0:0.02:0.08);
ytick(0:0.02:0.08);
xtickangle(0)
text(0.0025,0.075,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% clear col
% colormap(cmocean('ice',x_num))
% cbar=colorbar;
% caxis([x_nond(1),x_nond(end)])
% cbar.Label.FontSize = 25;
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$\mathrm{outer} \leftarrow \frac{x}{x_{b}} \rightarrow \mathrm{inner}$';
ax=gca;
%ax.YTickMode = "manual";
grid on 
%set(ax.XLabel,'Position',[ax.XLabel.Position(1),ax.XLabel.Position(2),ax.XLabel.Position(3)])
%set(ax.XLabel,'Position',[0.5000,-0.0690,0])




%============================================== panel b
h2=subplot("Position",pos(2,:));
clear col
%col = cmocean('thermal',x_num);
col_single = [.3 .3 .3];
col = repmat(col_single,5,1);

for xind = 1:x_num
scatter(A.kym2(:,xind),A.wkym2(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.kym3(:,xind),A.wkym3(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.kym4(:,xind),A.wkym4(:,xind),markersz,col(xind,:),'filled')
hold on 
end 
ylabel('$\hat{k}_{y0} \, (m^{-1})$','interpreter','latex','fontsize',16);
xlabel('$k_{y0} \, (m^{-1})$','interpreter','latex');
%title(['$r^2=$',num2str(rsquare_ky0,2)],'Interpreter','latex')
hold off 
niceplot_nobold(fig_fontsz)
grid off 
%axis equal
xlim([0 0.15])
ylim([0 0.15])
text(0.0045,0.14,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% clear col
% colormap(cmocean('thermal',x_num))
% cbar = colorbar;
% caxis([x_nond(1),x_nond(end)])
% cbar.Label.FontSize = 16;
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$x / x_{\mathrm{b}}$';
xtick(0:0.03:0.25)
ytick(0:0.03:0.25)
xtickangle(0)
plot_oneone 
grid on 
ax=gca;

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_one2one_G0_ky0.pdf'

%% AGU 2024 version
% plot_vfrc_scatter_col_slp(A.iS2s,A.wfit2.intw,1)
% plot_vfrc_scatter_col_slp(A.iS3s,A.wfit3.intw,2)
% plot_vfrc_scatter_col_slp(A.iS4s,A.wfit4.intw,3)
% ylabel('${\int_{0.001}^{0.2} \, S_{WB} \, dk_{y}} (s^{-4})$','interpreter','latex')
% xlabel('${\int_{0.001}^{0.2} \, S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, dk_{y}} (s^{-4})$','interpreter','latex')
% %title('With Cutoff','interpreter','latex')
% niceplot_nobold_nomintick(fig_fontsz)
% grid off
% %axis equal
% xlim([0 0.003])
% ylim([0 0.003])
% xtick(0:0.0005:0.003);
% ytick(0:0.0005:0.003);
% xtickangle(0)
% ax=gca;
% set(ax.XLabel,'Position',[0.00125,-5.0192e-04,-1])
% %ax.YTickMode = "manual";
% text(0.0001,0.0028,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% plot_oneone 
% grid on 
% %set(ax.XLabel,'Position',[ax.XLabel.Position(1),ax.XLabel.Position(2),ax.XLabel.Position(3)])
% %set(ax.XLabel,'Position',[0.5000,-0.0690,0])
% 
% 
% 
% 
% h2=subplot("Position",pos(2,:));
% plot_vfrc_scatter_col_slp(A.kym2,A.wkym2,1)
% plot_vfrc_scatter_col_slp(A.kym3,A.wkym3,2)
% plot_vfrc_scatter_col_slp(A.kym4,A.wkym4,3)
% ylabel('$\mathrm{Fit} \; k_{y0}  (\mathrm{cpm})$','interpreter','latex','fontsize',16);
% xlabel('$\mathrm{Model} \; k_{y0}  (\mathrm{cpm})$','interpreter','latex');
% %title('$k_{y0}$ from Data VS Weibull Fit','interpreter','latex')
% niceplot_nobold_nomintick(fig_fontsz)
% grid off 
% %axis equal
% xlim([0 0.15])
% ylim([0 0.15])
% xtick(0:0.03:0.15)
% ytick(0:0.03:0.15)
% xtickangle(0)
% plot_oneone 
% text(0.005,0.14,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% grid on 
% ax=gca;
% set(ax.XLabel,'Position',[0.075, -0.0261, -1])
% 
% 
% 
% 
% 
% set(gcf,'PaperUnit','centimeters')
% set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])
% 
% print -dpdf '/data1/bliu/vortpaper/fig_one2one_G0_ky0.pdf'

%close all
