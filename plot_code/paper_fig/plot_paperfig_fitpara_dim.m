% Bingchen Liu Mar 5,2025
% This code plot the nond G0 and nond ky0 scaling vs directional spread 




clear
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');
ky0=load('/data1/bliu/data/plotready_ky0_nond_10loc_2025.mat');



xsize=8;ysize=13.5;
x0=0.233; 
y0=0.1;  
dy=0.05;
xw=0.75; 
yw=0.395;
y1= y0+yw+dy; 
pos = [x0 y0 xw yw; x0 y1 xw yw];

axfont_sz = 15;
fig_fontsz= 10;
subfiglabel_fontsz = 15;
subfiglabel = {'(a)','(b)'};



markersz= 20;
x_num = length(G0.G0_nond.slp2(1,:));
clear col
x_nond = linspace(-0.75,-0.25,5);

loc_num=length(G0.G0_nond.slp2(1,:));
ds_10loc.slp2 =repmat(G0.ds.sigtb2,[loc_num,1])';
ds_10loc.slp3 =repmat(G0.ds.sigtb3,[loc_num,1])';
ds_10loc.slp4 =repmat(G0.ds.sigtb4,[loc_num,1])';


col = cmocean('thermal',x_num);

%%
figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:));
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0.G0_data.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o')
scatter(ds_10loc.slp3(:,xind),G0.G0_data.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square')
scatter(ds_10loc.slp4(:,xind),G0.G0_data.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^')
hold on 
end 
%xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
%title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
grid on 
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$\hat{G}_{0} \; (\mathrm{s}^{-4})$','interpreter','latex','fontsize',axfont_sz);
text(0.25,0.072,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
ax= gca;
set(gca,'XTickLabel',[])
%ax.YLim = [0 inf];
%xlim([0,15])
%ylim([0, 15])

hold off




col = cmocean('thermal',x_num);
subplot("Position",pos(1,:));
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0.ky0_data.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o')
scatter(ds_10loc.slp3(:,xind),ky0.ky0_data.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square')
scatter(ds_10loc.slp4(:,xind),ky0.ky0_data.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^')
hold on 
end 
grid on 
%title('Weibull Fit $k_{y0}$','interpreter','latex')
%xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$$\hat{k}_{y0} \; (\mathrm{m}^{-1})$$','interpreter','latex','fontsize',axfont_sz);
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex','fontsize',axfont_sz);
text(0.25,.11,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
hold off
clear col
colormap(cmocean('thermal',x_num))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.Label.FontSize = 12;
cbar.Label.Interpreter = 'latex';
cbar.Label.FontWeight = 'bold';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Ruler.TickLabelRotation=0;
cbar.Position = [0.379, 0.92,0.5,0.02];
cbar.AxisLocation='in';

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_dimG0ky0_ds.pdf'

close all