% Bingchen Liu Nov 29,2024
% This code plot the nond G0 and nond ky0 scaling vs directional spread 




clear
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');
ky0=load('/data1/bliu/data/plotready_ky0_nond_10loc_2025.mat');



xsize=8.5;ysize=13.5;
x0=0.22; 
y0=0.1;  
dy=0.05;
xw=0.76; 
yw=0.395;
y1= y0+yw+dy; 
pos = [x0 y0 xw yw; x0 y1 xw yw];

axfont_sz = 15;
fig_fontsz= 10;
subfiglabel_fontsz = 15;
font_rsq = 10;
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


%% double check r^2 using correlation calculation (not from fit function)

rsquare_G0 = G0.rsquare_G0(1,2)
rsquare_ky0 = ky0.rsquare_ky0(1,2)

['G0=',num2str(G0.f.A),'*dirspr+',num2str(G0.f.B)]
['ky0=',num2str(ky0.f.A),'*dirspr+',num2str(ky0.f.B)]

%%
figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:));
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0.G0_nond.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),G0.G0_nond.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),G0.G0_nond.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
end 
%xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
%title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
plot(G0.ds_model,G0.G0_model,'LineWidth',3,'Color','r','LineStyle','--')
grid on 
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$\tilde{G}_{0}=\frac{\hat{G}_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \beta$','interpreter','latex','fontsize',axfont_sz);
text(0.25,11,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
%text(11,0.1,['$r^2 = $',num2str(G0.gof.rsquare,2)],'FontSize',font_rsq,'Interpreter','latex')
ax= gca;
set(gca,'XTickLabel',[])
%ax.YLim = [0 inf];
xlim([0,15])
%ylim([0, 15])

hold off




col = cmocean('thermal',x_num);
subplot("Position",pos(1,:));
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0.ky0_nond.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),ky0.ky0_nond.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),ky0.ky0_nond.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
end 
plot(ky0.ds_model,ky0.ky0_model,'LineWidth',3,'Color','r','LineStyle','--')
grid on 
%title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$$\tilde{k}_{y0} = \hat{k}_{y0} \, h \, \mathrm{Ir}_{b}$$','interpreter','latex','fontsize',axfont_sz);
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex','fontsize',axfont_sz);
text(0.25,.023,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
%text(11,0.00025,['$r^2 = $',num2str(ky0.gof.rsquare,2)],'FontSize',font_rsq,'Interpreter','latex')
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
cbar.Position = [0.3, 0.92,0.5,0.02];
cbar.AxisLocation='in';
ax=gca;
ax.YAxis.Exponent = -2;
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_nondG0ky0_ds.pdf'

close all