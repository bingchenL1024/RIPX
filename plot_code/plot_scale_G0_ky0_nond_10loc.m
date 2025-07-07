% Bingchen Liu Nov 29,2024
% This code plot the nond G0 and nond ky0 scaling vs directional spread 
% modified from 'plot_scal_nond' --> without saving actual figure 



clearvars -except expo


G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');
ky0=load('/data1/bliu/data/plotready_ky0_nond_10loc_2025.mat');

loc_num=length(G0.G0_nond.slp2(1,:));

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


markersz= 80;
x_num = length(G0.G0_nond.slp2(1,:));
clear col
x_nond = linspace(-0.75,-0.25,5);

ds_10loc.slp2 =repmat(G0.ds.sigtb2,[loc_num,1])';
ds_10loc.slp3 =repmat(G0.ds.sigtb3,[loc_num,1])';
ds_10loc.slp4 =repmat(G0.ds.sigtb4,[loc_num,1])';


col = cmocean('ice',x_num);

figure()
subplot(211);
plot(G0.ds_model,G0.G0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0.G0_nond.slp2(:,xind),markersz,col(xind,:),'filled')
scatter(ds_10loc.slp3(:,xind),G0.G0_nond.slp3(:,xind),markersz,col(xind,:),'filled')
scatter(ds_10loc.slp4(:,xind),G0.G0_nond.slp4(:,xind),markersz,col(xind,:),'filled')
end 
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
grid on 
clear col
colormap(cmocean('ice',x_num))
cbar=colorbar;
caxis([x_nond(1),x_nond(end)])
cbar.Label.FontSize = 25;
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\mathrm{outer} \leftarrow \frac{x}{x_{b}} \rightarrow \mathrm{inner}$';
hold off

col = cmocean('ice',x_num);

subplot(212);
plot(ky0.ds_model,ky0.ky0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0.ky0_nond.slp2(:,xind),markersz,col(xind,:),'filled')
scatter(ds_10loc.slp3(:,xind),ky0.ky0_nond.slp3(:,xind),markersz,col(xind,:),'filled')
scatter(ds_10loc.slp4(:,xind),ky0.ky0_nond.slp4(:,xind),markersz,col(xind,:),'filled')
end 
%title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(fig_fontsz)
title(['r^2 =',num2str(ky0.gof.rsquare)])
ylabel('$$k_{y0} \, h \, \beta$$','interpreter','latex','fontsize',axfont_sz);
xlabel('$\sigma_\theta$~(deg)','interpreter','latex','fontsize',axfont_sz);
clear col
colormap(cmocean('ice',x_num))
cbar=colorbar;
caxis([x_nond(1),x_nond(end)])
cbar.Label.FontSize = 25;
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\mathrm{outer} \leftarrow \frac{x}{x_{b}} \rightarrow \mathrm{inner}$';

grid on 
hold off


