% Bingchen Liu June 4, 2025
% This code plot the figure for slides for talk (reduced complexicity)
% file saved to /bliu/analysis_code/talk_fig


%% 2D vort and forcing field 
clear
close all
load('/data1/bliu/data/snap_vort_curlF_field_run14.mat')
load('/data1/bliu/data/SS_raw.mat')
ind=14;


xsize=12;ysize=13;
x0=0.162; 
y0=0.102;  
dy=0.03;
xw=0.75; 
yw=xw*(xsize/ysize)*(210/350);
y1= y0+yw+dy; 
pos = [x0 y0 xw yw; x0 y1 xw yw];


xb = S(ind).xb;
% vort_snap = flip(vort_snap,2);
% curlF_snap = flip(curlF_snap,2);

dim= size(curlF_snap');
% x = 0:dim(1)-1;
% y = 0:dim(2)-1; 

x = -(dim(1)-1):0; %modified
y = 0:dim(2)-1; %modified
[x_grid,y_grid] = meshgrid(x,y);
% for t = 2000:2100;



fig_fontsz= 10;
subfiglabel_fontsz = 16;
axfont_sz = 14;
col_fontsz= 14;
subfiglabel = {'(a)','(b)'};

% ======================= plotting 

figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:));
h1= pcolorcen(y_grid,x_grid,vort_snap);
yline(-xb,'k--','LineWidth',2)
%xline(10,'k-','LineWidth',2)
col=colorbar;
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-0.12,0.12])
col.Label.Interpreter = 'latex';
col.Label.String = '$\omega \, (s^{-1})$';
col.Label.FontSize = col_fontsz;
col.Label.FontWeight = 'bold';
col.Ticks = linspace(-0.1,0.1,5);
%xlabel('$y$ (m)','Interpreter','latex')
%title('$\nabla \times \vec{u}$','Interpreter','latex')
%title([runpara.wave,runpara.bath])
% width= 22;
% height = 18;
% set(gcf,'Units','inches','Position',[0,0,width,height]);
% set(gcf,'visible','off');
%ylim([dim(1)-100,dim(1)])
%axis equal
% ylim([0,dim(1)])
% xlim([0,350])
% xtick(0:50:350)
% ytick(0:30:90)

 ylim([-dim(1),0]) %modified
 xlim([0,350]) %modified
 ytick(-250:50:0) %modified
 set(gca,'YDir','reverse') %modified

%ax = gca;
niceplot_nobold(fig_fontsz)
text(15,-240,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
%text(150,-5,'Beach','FontSize',subfiglabel_fontsz)

ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
% ax.XTick = [];
%ax.XTickLabel = [];
xtick(0:50:300) %modified
set(gca,'XTickLabel',[])


subplot("Position",pos(1,:));
pcolorcen(y_grid,x_grid,curlF_snap);
yline(-xb,'k--','LineWidth',2)
%xline(10,'k-','LineWidth',2)
col=colorbar;
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-0.3,0.3])
%col.Limits = [-0.35,0.35];
col.Label.Interpreter = 'latex';
col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = col_fontsz;
%col.Label.FontWeight = 'bold';
col.Ticks = linspace(-0.3,0.3,5);
%col.TickLabels=-0.5:0.1:0.5;
%title('$\nabla \times \vec{u}$','Interpreter','latex')
% width= 25;
% height = 10;
% set(gcf,'Units','inches','Position',[0,0,width,height]);
% set(gcf,'visible','off');
%set(gca,'YDir','reverse')
%axis equal
ylim([-dim(1),0]) %modified
xlim([0,350]) %modified
xtick(0:50:300) %modified
ytick(-250:50:0) %modified
set(gca,'YDir','reverse') %modified
niceplot_nobold(fig_fontsz)
text(15,-240,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
%text(150,-5,'Beach','FontSize',subfiglabel_fontsz)
ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
xlabel('$y$ (m)','Interpreter','latex','FontSize',axfont_sz)
ax2=gca;
set(ax2.XAxis,'TickDirection','out')


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/analysis_code/talk_fig/fig_snap2d_vortcurlF_field.pdf'
close all 

%%



clear
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');
ky0=load('/data1/bliu/data/plotready_ky0_nond_10loc_2025.mat');



xsize=8;ysize=13.5;
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


figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:));

plot(G0.ds_model,G0.G0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0.G0_nond.slp2(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),G0.G0_nond.slp3(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),G0.G0_nond.slp4(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
end 
%xlabel('$\sigma_\theta$~(deg)','interpreter','latex');
%title('$G_{0} = \sqrt{\int^{k_{max}} \, S_{\mathrm{WB}} \, dk_{y}}$','interpreter','latex')
grid on 
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$\frac{\hat{G}_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \, \mathrm{Ir}_{\mathrm{b}}$','interpreter','latex','fontsize',axfont_sz);
text(0.25,11,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
ax= gca;
set(gca,'XTickLabel',[])
%ax.YLim = [0 inf];
xlim([0,15])
%ylim([0, 15])

hold off




col = cmocean('thermal',x_num);
subplot("Position",pos(1,:));
plot(ky0.ds_model,ky0.ky0_model,'LineWidth',3,'Color','k','LineStyle','--')
hold on 
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0.ky0_nond.slp2(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),ky0.ky0_nond.slp3(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),ky0.ky0_nond.slp4(:,xind),markersz,col(xind,:),'filled','MarkerEdgeColor','k')
end 
grid on 
%title('Weibull Fit $k_{y0}$','interpreter','latex')
xlim([0,15])
%ylim([0,0.1])
niceplot_nobold_nomintick(fig_fontsz)
ylabel('$$\hat{k}_{y0} \, h \, \mathrm{Ir}_{\mathrm{b}}$$','interpreter','latex','fontsize',axfont_sz);
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex','fontsize',axfont_sz);
text(0.25,.023,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
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
ax=gca;
ax.YAxis.Exponent = -2;
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/analysis_code/talk_fig/fig_nondG0ky0_ds.pdf'

close all
