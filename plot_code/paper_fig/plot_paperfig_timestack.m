% Bingchen Liu Feb 5, 2025
% This code plot the time stack of vort and curlF for paper fig
clear
load('/data1/bliu/data/tstack_curlF_run14.mat') %'get_timestack'
load('/data1/bliu/data/tstack_vort_run14.mat') %'get_timestack'
load('/data1/bliu/data/SS_raw.mat')


label_fontsz = 16;
col_fontsz = 15;
subfiglabel = {'(a)','(b)'};
subfiglabel_fontsz = 15;
fig_fontsz = 12;

xsize=12;ysize=12;
x0=0.16; 
y0=0.115;  
dy=0.03;
xw=0.742; 
yw=xw*(xsize/ysize)*(195/350);
y1= y0+yw+dy; 
pos = [x0 y0 xw yw; x0 y1 xw yw];

%%
ind=14;
xb = S(ind).xb;

dim = size(curlF_tstack);

x = -(dim(1)-1):0;
t = 0:dim(2)-1; %in seconds
[t_grid,x_grid] = meshgrid(t,x);


%%
figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(2,:));
pcolorcen(t_grid,x_grid,vort_tstack);
hold on 
yline(-xb,'k--','LineWidth',2)
hold off


col=colorbar;
cmocean('balance');
caxis([-0.12,0.12])
set(gca,'YDir','reverse') %modified
col.Label.Interpreter = 'latex';
col.Label.String = '$\omega \, (s^{-1})$';
col.Label.FontSize = col_fontsz;
col.Label.FontWeight = 'bold';
col.Ticks = linspace(-0.1,0.1,5);

ylabel('$x$ (m)','Interpreter','latex','FontSize',label_fontsz)
%xlabel('$t$ (s)','Interpreter','latex','FontSize',label_fontsz)
xtick(0:25:75) %modified
xlim([0,80])
set(gca,'XTickLabel',[])

niceplot_nobold(fig_fontsz)
text(1,-240,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')


%%
subplot("Position",pos(1,:));
pcolorcen(t_grid,x_grid,curlF_tstack);
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

ylabel('$x$ (m)','Interpreter','latex','FontSize',label_fontsz)
xlabel('$t$ (s)','Interpreter','latex','FontSize',label_fontsz)
xtick(0:25:75) %modified
xlim([0,80])

niceplot_nobold(fig_fontsz)
text(1,-240,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/vortpaper/fig_Hovmoller.pdf'
close all 
