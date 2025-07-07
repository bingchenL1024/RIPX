% Bingchen Liu Jan 23, 2025
% This plot code plot one-to-one plot for G0 and ky0 compare model data
% with Weibull fit 


load('/data1/bliu/data/plotready_G0_nond_10loc_2025')
load('/data1/bliu/data/G0_nond_10loc_allsk_2025')


fig_fontsize = 22;
markersz= 80;
x_num = length(A.iS2s(1,:));
clear col
col = cmocean('ice',x_num);
x_nond = -1:0.1:-0.1;

figure()
subplot(121)
for xind = 1:x_num
scatter(A.iS2s(:,xind),A.wfit2.intw(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.iS3s(:,xind),A.wfit3.intw(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.iS4s(:,xind),A.wfit4.intw(:,xind),markersz,col(xind,:),'filled')
hold on 
end 
plot_oneone 
niceplot_nobold(fig_fontsize)
hold off 
ylabel('${\int_{0.001}^{0.2} \, S_{WB} \, dk_{y}} (s^{-4})$','interpreter','latex')
xlabel('${\int_{0.001}^{0.2} \, S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, dk_{y}} (s^{-4})$','interpreter','latex')
%axis equal
%xlim([0 0.003])
%ylim([0 0.003])
%xtick(0:0.0005:0.003);
%ytick(0:0.0005:0.003);
xtickangle(0)
clear col
colormap(cmocean('ice',x_num))
cbar=colorbar;
caxis([x_nond(1),x_nond(end)])
cbar.Label.FontSize = 25;
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\mathrm{outer} \leftarrow \frac{x}{x_{b}} \rightarrow \mathrm{inner}$';
ax=gca;
%ax.YTickMode = "manual";
grid on 
%set(ax.XLabel,'Position',[ax.XLabel.Position(1),ax.XLabel.Position(2),ax.XLabel.Position(3)])
%set(ax.XLabel,'Position',[0.5000,-0.0690,0])
title(['Expo = ',num2str(expo)])



clear col
col = cmocean('ice',x_num);

subplot(122)
for xind = 1:x_num
scatter(A.kym2(:,xind),A.wkym2(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.kym3(:,xind),A.wkym3(:,xind),markersz,col(xind,:),'filled')
hold on 
scatter(A.kym4(:,xind),A.wkym4(:,xind),markersz,col(xind,:),'filled')
hold on 
ylabel('$\mathrm{Fit} \; k_{y0}  (\mathrm{cpm})$','interpreter','latex','fontsize',16);
xlabel('$\mathrm{Model} \; k_{y0}  (\mathrm{cpm})$','interpreter','latex');
%title('$k_{y0}$ from Data VS Weibull Fit','interpreter','latex')
end 
hold off 
niceplot_nobold(fig_fontsize)
grid off 
%axis equal
%xlim([0 0.15])
%ylim([0 0.15])
clear col
colormap(cmocean('ice',x_num))
cbar = colorbar;
caxis([x_nond(1),x_nond(end)])
cbar.Label.FontSize = 25;
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\mathrm{outer} \leftarrow \frac{x}{x_{b}} \rightarrow \mathrm{inner}$';
niceplot_nobold_nomintick(18);
xtick(0:0.03:0.25)
ytick(0:0.03:0.25)
xtickangle(0)
plot_oneone 
grid on 
ax=gca;

