% Bingchen Liu Nov 27, 2024
% This code plot a snopshot of vorticity and curlF field for paper


%% negative cross-shore coordinate  

clear
close all
load('/data1/bliu/data/snap_vort_curlF_field_run14.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/tstack_curlF_run14.mat') %'get_timestack'



%%
xsize=12;ysize=20;
x0=0.162; 
y0=0.0635;  
dy=0.03;
xw=0.75; 
yw=xw*(xsize/ysize)*(205/350); %ratio of xlim and ylim
y1= y0+yw+dy+y0; y2 = y1+dy+yw;
pos = [x0 y0 xw yw; x0 y1 xw yw;x0 y2 xw yw];

%% datat info
ind=14;

dim_Hov = size(curlF_tstack);

x_Hov = -(dim_Hov(1)-1):0;
t_Hov = 0:dim_Hov(2)-1; %in seconds
[t_grid_Hov,x_grid_Hov] = meshgrid(t_Hov,x_Hov);
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

box_pos= [85,-170,45,20] ;

fig_fontsz= 10;
subfiglabel_fontsz = 16;
axfont_sz = 14;
col_fontsz= 14;
subfiglabel = {'(a)','(b)','(c)'};

% ======================= plotting 

figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


subplot("Position",pos(3,:));
h1= pcolorcen(y_grid,x_grid,vort_snap);
yline(-xb,'k--','LineWidth',2)
%xline(10,'k-.','LineWidth',2)
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


subplot("Position",pos(2,:));
pcolorcen(y_grid,x_grid,curlF_snap);
yline(-xb,'k--','LineWidth',2)
xline(10,'k-.','LineWidth',2)
hold on 
rectangle('Position',box_pos,'EdgeColor','k','LineStyle','--','LineWidth',2)
hold off 
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

subplot("Position",pos(1,:));
pcolorcen(t_grid_Hov,x_grid_Hov,curlF_tstack);
hold on 
yline(-xb,'k--','LineWidth',2)
hold off
col=colorbar;
cmocean('balance');
caxis([-0.3,0.3])
set(gca,'YDir','reverse') %modified
col.Label.Interpreter = 'latex';
col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = col_fontsz;
col.Label.FontWeight = 'bold';
col.Ticks = linspace(-0.3,0.3,5);
niceplot_nobold(fig_fontsz)
ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
xlabel('$t$ (s)','Interpreter','latex','FontSize',axfont_sz)
xtick(0:25:75) %modified
xlim([0,80])
text(3.4,-240,subfiglabel{3},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/vortpaper/fig_snap2d_vortcurlF_field.pdf'
close all 

%% positive cross-shore coordinate 
% clear
% close all
% load('/data1/bliu/data/snap_vort_curlF_field')
% load('/data1/bliu/data/SS_raw.mat')
% 
% xsize=14;ysize=8;
% x0=0.12; 
% y0=0.17;  
% dy=0.05;
% xw=0.8; 
% yw=xw*(xsize/ysize)*(94/350);
% y1= y0+yw+dy; 
% pos = [x0 y0 xw yw; x0 y1 xw yw];
% 
% 
% ind=44;
% xb = S(ind).xb;
% vort_snap = flip(vort_snap,2);
% curlF_snap = flip(curlF_snap,2);
% 
% dim= size(curlF_snap');
% x = 0:dim(1)-1;
% y = 0:dim(2)-1; 
% 
% 
% [x_grid,y_grid] = meshgrid(x,y);
% % for t = 2000:2100;
% 
% 
% 
% fig_fontsz= 10;
% subfiglabel_fontsz = 14;
% axfont_sz = 14;
% col_fontsz= 14;
% subfiglabel = {'(a)','(b)'};
% 
% % ======================= plotting 
% 
% figure()
% 
% set(gcf,'PaperUnit','centimeters')
% set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])
% 
% 
% subplot("Position",pos(2,:));
% h1= pcolorcen(y_grid,x_grid,vort_snap);
% yline(xb,'k--','LineWidth',2)
% col=colorbar;
% cmocean('balance');
% %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
% caxis([-0.1,0.1])
% col.Label.Interpreter = 'latex';
% col.Label.String = '$\omega \, (s^{-1})$';
% col.Label.FontSize = col_fontsz;
% col.Label.FontWeight = 'bold';
% col.Ticks = linspace(-0.1,0.1,5);
% %xlabel('$y$ (m)','Interpreter','latex')
% %title('$\nabla \times \vec{u}$','Interpreter','latex')
% %title([runpara.wave,runpara.bath])
% % width= 22;
% % height = 18;
% % set(gcf,'Units','inches','Position',[0,0,width,height]);
% % set(gcf,'visible','off');
% %ylim([dim(1)-100,dim(1)])
% %axis equal
% ylim([0,dim(1)])
% xlim([0,350])
% xtick(0:50:350)
% ytick(0:30:90)
% 
%  % ylim([-dim(1),0]) %modified
%  % xlim([0,350]) %modified
%  % xtick(0:50:350) %modified
%  % ytick(-90:30:0) %modified
%  % set(gca,'YDir','reverse') %modified
% 
% %ax = gca;
% niceplot_nobold(fig_fontsz)
% text(1,86,subfiglabel{1},'FontSize',subfiglabel_fontsz)
% ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
% % ax.XTick = [];
% %ax.XTickLabel = [];
% set(gca,'XTickLabel',[])
% 
% 
% subplot("Position",pos(1,:));
% pcolorcen(y_grid,x_grid,curlF_snap);
% yline(xb,'k--','LineWidth',2)
% col=colorbar;
% cmocean('balance');
% %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
% caxis([-0.5,0.5])
% %col.Limits = [-0.35,0.35];
% col.Label.Interpreter = 'latex';
% col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
% col.Label.FontSize = col_fontsz;
% %col.Label.FontWeight = 'bold';
% col.Ticks = linspace(-0.5,0.5,5);
% %col.TickLabels=-0.5:0.1:0.5;
% %title('$\nabla \times \vec{u}$','Interpreter','latex')
% % width= 25;
% % height = 10;
% % set(gcf,'Units','inches','Position',[0,0,width,height]);
% % set(gcf,'visible','off');
% %set(gca,'YDir','reverse')
% %axis equal
% ylim([0,dim(1)])
% xlim([0,350])
% xtick(0:50:340)
% ytick(0:30:90)
% niceplot_nobold(fig_fontsz)
% text(1,86,subfiglabel{2},'FontSize',subfiglabel_fontsz)
% ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
% xlabel('$y$ (m)','Interpreter','latex','FontSize',axfont_sz)
% ax2=gca;
% set(ax2.XAxis,'TickDirection','out')
% 
% 
% set(gcf,'PaperUnit','centimeters')
% set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])
% 
% % print -dpdf -painters '/data1/bliu/vortpaper/fig_snap2d_vortcurlF_field.pdf'
% % close all 