% Bingchen Liu July 5, 2025
% This code plot the 2 panels plot (snap + timestack) using the
% parameterized forcing 


clear
close all
%load('/data1/bliu/data/parameterization_example')
load('/data1/bliu/data/parameterization_example_realistic.mat')
load('/data1/bliu/data/SS_raw.mat')

X= X_full;
Y=Y_full;
%%
xsize=12;ysize=18;
x0=0.1536; 
y0=0.075;  
dy=0.03;
xw=0.75; 
yw=xw*(xsize/ysize)*0.8; %ratio of xlim and ylim
y1= y0+yw+dy+y0; y2 = y1+dy+yw;
pos = [x0 y0 xw yw; x0 y1 xw yw];



%% datat info
ind=14;
xb = S(ind).xb;

dim_Hov = size(curlF_tstack);

x_Hov = -(dim_Hov(1)-1):0;
t_Hov = 0:dim_Hov(2)-1; %in seconds
[t_grid_Hov,x_grid_Hov] = meshgrid(t_Hov,x_Hov);
% vort_snap = flip(vort_snap,2);
% curlF_snap = flip(curlF_snap,2);

dim= size(curlF_snap');
% x = 0:dim(1)-1;
% y = 0:dim(2)-1; 

x = -(dim(1)-1):0; %modified
y = 0:dim(2)-1; %modified
[x_grid,y_grid] = meshgrid(x,y);
% for t = 2000:2100;

clim = 0.3; %colorbar limit for cFbr



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
pcolorcen(y_grid,x_grid,curlF_snap);
xline(10,'k-.','LineWidth',2)
yline(-xb,'k--','LineWidth',2)
hold off 
col=colorbar;
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-clim,clim])
%col.Limits = [-0.35,0.35];
col.Label.Interpreter = 'latex';
col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = col_fontsz;
%col.Label.FontWeight = 'bold';
col.Ticks = linspace(-clim,clim,5);
ylim([-dim(1),0]) %modified
xlim([0,300]) %modified
xtick(0:100:300) %modified
ytick(-300:100:0) %modified
set(gca,'YDir','reverse') %modified
niceplot_nobold(fig_fontsz)
text(10,-280,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
xlabel('$y$ (m)','Interpreter','latex','FontSize',axfont_sz)
ax2=gca;
set(ax2.XAxis,'TickDirection','out')
set(gca,'Layer','top')

subplot("Position",pos(1,:));
pcolorcen(t_grid_Hov,x_grid_Hov,curlF_tstack);
hold on 
yline(-xb,'k--','LineWidth',2)
hold off
col=colorbar;
cmocean('balance');
caxis([-clim,clim])
set(gca,'YDir','reverse') %modified
col.Label.Interpreter = 'latex';
col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = col_fontsz;
col.Label.FontWeight = 'bold';
col.Ticks = linspace(-clim,clim,5);
niceplot_nobold(fig_fontsz)
ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
xlabel('$t$ (s)','Interpreter','latex','FontSize',axfont_sz)
xtick(0:10:70) %modified
ytick(-300:100:0) %modified
ylim([-dim(1),0]) %modified
xlim([0,70])
text(3,-280,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
set(gca,'Layer','top')

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/vortpaper/fig_exampleF_2panels.pdf'