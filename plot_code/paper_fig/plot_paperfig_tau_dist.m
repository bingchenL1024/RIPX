% Bingchen Liu April 2, 2025
% This code plot the distribution of \tau for paper fig

load('/data1/bliu/data/cxt_fitanalysis_max_dxwidth.mat')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'

%% 1 panel 
upperlim = 1.75;
lowlim = 0.25;
bin_num=10;
bin_edge = linspace(lowlim,upperlim,bin_num+1);
tick = 0.25:0.25:upperlim;

xsize=8.5;ysize=8;
x0=0.15; 
y0=0.15;  
dx=0.12;
xw=0.75; yw=xw*(xsize/ysize);
x1= x0+xw+dx; 
pos = [x0 y0 xw yw];

axfont_sz = 15;
fig_fontsz= 10;
mksz = 66;
mkwidth = 1.2;
subfiglabel_fontsz = 14;

fig=figure();
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

h1=histogram(a_tot_all,'BinEdges',bin_edge,'Normalization','probability','FaceAlpha',1);
hold on 
scatter(t_deco_1mres_stats.mean,0,mksz,'d','MarkerEdgeColor','m','LineWidth',mkwidth)
hold on 

xlabel('$\tau \, (s)$','Interpreter','latex')
%ylabel('Fraction of occurences')
xticks(tick)
scatter(t_deco_1mres_stats.med,0,mksz,'^','MarkerEdgeColor','m','LineWidth',mkwidth)
set(gca,'Layer','bottom')
xtickangle(0)
hold off
grid on 
ylabel('Fraction of occurences')
niceplot_nobold(fig_fontsz)


print -dpdf '/data1/bliu/vortpaper/fig_tau_dist.pdf'

close all





%% 2 panels 
% bin_num=10;
% bin_edge = linspace(0.35,1.5,bin_num+1);
% tick = 0.5:0.25:1.5;
% 
% xsize=7.5;ysize=13.5;
% x0=0.22; 
% y0=0.1;  
% dy=0.05;
% xw=0.76; 
% yw=0.395;
% y1= y0+yw+dy; 
% pos = [x0 y0 xw yw; x0 y1 xw yw];
% 
% axfont_sz = 15;
% fig_fontsz= 10;
% mksz = 66;
% mkwidth = 1.2;
% subfiglabel_fontsz = 14;
% subfiglabel = {'(a)','(b)'};
% 
% fig=figure();
% set(gcf,'PaperUnit','centimeters')
% set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])
% 
% subplot("Position",pos(2,:));
% h1=histogram(a_tot_all,'BinEdges',bin_edge,'Normalization','probability','FaceAlpha',1);
% hold on 
% scatter(t_deco_1mres_stats.mean,0,mksz,'d','MarkerEdgeColor','k','LineWidth',mkwidth)
% hold on 
% scatter(t_deco_1mres_stats.med,0,mksz,'^','MarkerEdgeColor','k','LineWidth',mkwidth)
% 
% %xlabel('$\tau$','Interpreter','latex')
% %ylabel('Fraction of occurences')
% xticks(tick)
% grid on 
% set(gca,'XTickLabel',[])
% niceplot_nobold(fig_fontsz)
% text(0.32,0.234,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% 
% 
% 
% 
% subplot("Position",pos(1,:));
% h2=histogram(t_decoscal_tot_5loc,'BinEdges',bin_edge,'Normalization','probability','FaceAlpha',1);
% hold on 
% scatter(fitpara_5loc.t_mean,0,mksz,'d','MarkerEdgeColor','k','LineWidth',mkwidth)
% hold on 
% scatter(fitpara_5loc.t_med,0,mksz,'^','MarkerEdgeColor','k','LineWidth',mkwidth)
% 
% xlabel('$\tau (s)$','Interpreter','latex')
% %ylabel('Fraction of occurences')
% xticks(tick)
% niceplot_nobold(fig_fontsz)
% text(0.32,0.283,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
% grid on 
% 
% han = axes(fig,'Visible','off');
% han.YLabel.Visible= 'on';
% ylabel(han,'Fraction of occurences')

