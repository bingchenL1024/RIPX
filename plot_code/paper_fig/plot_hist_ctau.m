% Bingchen Liu July 1, 2025 
% This code plots the historgram of c*\tau 


load('/data1/bliu/data/cxt_fitanalysis_max_dxwidth.mat') %'get_cxt_fitanalysis_max_dxwidth'
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'

%% 1 panel 
upperlim = 8;
lowlim = 0;
bin_num=16;
bin_edge = linspace(lowlim,upperlim,bin_num+1);
tick = lowlim:1:upperlim;

xsize=9;ysize=7;
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

h1=histogram(ctau_all,'BinEdges',bin_edge,'Normalization','probability','FaceAlpha',1);
hold on 
scatter(ctau_stats.mean,0,mksz,'d','MarkerEdgeColor','m','LineWidth',mkwidth)
hold on 
scatter(ctau_stats.med,0,mksz,'^','MarkerEdgeColor','m','LineWidth',mkwidth)

xlabel('$\hat{c}\hat{\tau} \, (\mathrm{m})$','Interpreter','latex')
%ylabel('Fraction of occurences')
xticks(tick)
set(gca,'Layer','bottom')
xtickangle(0)
hold off
grid on 
ylabel('Fraction of occurences')
niceplot_nobold(fig_fontsz)


print -dpdf '/data1/bliu/vortpaper/fig_ctau_dist.pdf'
