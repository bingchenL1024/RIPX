% Bingchen Liu June 29,2025
% this code plots the blowed up forcing structure for paper discussion
% section


%% negative cross-shore coordinate  

clear
close all
load('/data1/bliu/data/snap_vort_curlF_field_run14.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/tstack_curlF_run14.mat') %'get_timestack'

%% fig dimension 
xsize=15;ysize=8;
x0=0.12; 
y0=0.2;  
dx=0.12;
xw=0.85; yw=xw*((ysize+4.5)/xsize);
x1= x0+xw+dx; 
pos = [x0 y0 xw yw];

cpos_x = 0.35; %colorbar cocation
cpos_y=0.77;

fig_fontsz= 10;
subfiglabel_fontsz = 15;
cbar_fontsize =10;
cbar_ticksz= 6;

markersz= 20;

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

ydomain_ind = 86:131;%x and y domain of the Zoom in region for the focing 
xdomain_ind = 89:109; 
xdomain = x_grid(1,xdomain_ind);
ydomain = y_grid(ydomain_ind,1);

fig_fontsz= 10;
subfiglabel_fontsz = 16;
axfont_sz = 14;
col_fontsz= 12;
subfiglabel = {'(a)','(b)','(c)'};

% ======================= plotting 

figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])



subplot("Position",pos(1,:));
pcolorcen(y_grid(ydomain_ind,xdomain_ind),x_grid(ydomain_ind,xdomain_ind),curlF_snap(ydomain_ind,xdomain_ind));
colormap(cmocean('balance'))
col=colorbar('south');

ylim([xdomain(1),xdomain(end)]) %modified
xlim([ydomain(1),ydomain(end)]) %modified
xtick(85:5:130) %modified
ytick(-170:5:-150) %modified
set(gca,'YDir','reverse') %modified
niceplot_nobold(fig_fontsz)
text(15,-240,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
%text(150,-5,'Beach','FontSize',subfiglabel_fontsz)
ylabel('$x$ (m)','Interpreter','latex','FontSize',axfont_sz)
xlabel('$y$ (m)','Interpreter','latex','FontSize',axfont_sz)
ax2=gca;
set(ax2.XAxis,'TickDirection','out')

%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-0.3,0.3])
%col.Limits = [-0.35,0.35];
col.Label.Interpreter = 'latex';
%col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = col_fontsz;
%col.Label.FontWeight = 'bold';
col.Ticks = linspace(-0.3,0.3,5);
col.Position = [cpos_x,cpos_y,0.39,0.05];
set(col,'XAxisLocation','top')
col.Ruler.TickLabelRotation=0;
col.Label.Rotation = 0;
%set(col,'XTick',[-1,-1,-,])

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/vortpaper/fig_forcing_detail.pdf'
close all 