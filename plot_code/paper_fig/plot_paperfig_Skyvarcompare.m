% Bingchen Liu Jan 30, 2025
% This code plot paper figure to compare total spectra variance to spectra
% variance with cutoff

clear
close all
% load('/data1/bliu/data/Sky_nond')
% load('/data1/bliu/data/Sky_WBfit_qced') 

load('data1/bliu/data/Sky_nond_10loc_1p4.mat') %'get_weibul_qced_10loc'
load('/data1/bliu/data/Sky_binmean_10loc_1p4')


xsize=8;ysize=8;
x0=0.25; 
y0=0.15;  
dx=0.12;
xw=0.69; yw=xw*(xsize/ysize);
x1= x0+xw+dx; 
pos = [x0 y0 xw yw];

fig_fontsz= 10;
subfiglabel_fontsz = 15;
subfiglabel = {'(a)','(b)'};

markersz= 20;

%% stats ('get_weibul_qced_10loc')
stats.rsquare_G0full
stats.rmse_G0full 

%%

x_num = length(A.iS2s(1,:));
x_nond = linspace(-0.75,-0.25,5);
%col = cmocean('thermal',x_num);
col_single = [.3 .3 .3];
col = repmat(col_single,5,1);

figure()
h1=subplot("Position",pos(1,:));
for xind = 1:x_num
scatter((A.iS2tot(:,xind)).^(0.5),(A.iS2s(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
scatter((A.iS3tot(:,xind)).^(0.5),(A.iS3s(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
scatter((A.iS4tot(:,xind)).^(0.5),(A.iS4s(:,xind)).^(0.5),markersz,col(xind,:),'filled')
hold on 
end 
plot_oneone 
niceplot_nobold(fig_fontsz)
hold off 
xlabel('$G_{\mathrm{full}} \, (s^{-2})$','interpreter','latex')
ylabel('$G_{0} \, (s^{-2})$','Interpreter','latex') %'${\int_{0.001}^{0.2} \, S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, dk_{y}} (s^{-4})$','interpreter','latex')
%title(['$r^2=$',num2str(rsquare_G0var,3)],'Interpreter','latex')
xlim([0 0.08])
ylim([0 0.08])
ax=gca;
%set(ax.XLabel,'Position',[0.00225,-6.5192e-04,-1])%xtick(0:0.0005:0.003);
ytick(0:0.02:0.08);
xtick(0:0.02:0.08);
xtickangle(0)
% clear col
% colormap(cmocean('thermal',x_num))
% cbar = colorbar;
% caxis([x_nond(1),x_nond(end)])
% cbar.Label.FontSize = 16;
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$x / x_{\mathrm{b}}$';
% xtick(0:0.03:0.25)
% ytick(0:0.03:0.25)
% xtickangle(0)
grid on 


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_Skyvarcompare.pdf'