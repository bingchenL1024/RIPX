% Bingchen Liu June 30, 2025
% This code plot the example forcing for the paper discussion section 
clear
load('/data1/bliu/data/parameterization_example')


%%
xsize=15;ysize=8;
x0=0.13; 
y0=0.16;  
dy=0.06;
dx=0.15;
xw=0.335; 
yw=0.345; %ratio of xlim and ylim
y1= y0+yw+dy; y2 = y1+dy+yw;
x1 = x0+xw+dx;
pos = [x0 y0 xw yw; x0 y1 xw yw;x1 y0 xw+0.05 2*yw+dy];


%
tind= 10;
xmax= -200;

fig_fontsz= 10;
subfiglabel_fontsz = 16;
axfont_sz = 12;
col_fontsz= 14;
linewd = 1.1;
subfiglabel = {'(a)','(b)','(c)'};
sublabelposx=-195;

%
close all
figure()
subplot("Position",pos(2,:));
plot(-X(1,:),G0(1,:,tind),'Color','k','LineWidth',linewd)
xlim([xmax,0])
ylim([0,0.08])
ytick(0:0.02:0.08)
set(gca,'XTickLabel',[])
ylabel('$G_{0} \, (\mathrm{m}^{-2})$','Interpreter','latex')
niceplot_nobold(axfont_sz)
text(sublabelposx,0.071,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

subplot("Position",pos(1,:));
plot(-X(1,:),W(1,:,tind),'Color','k','LineWidth',linewd)
xlabel('$x \, (\mathrm{m})$','Interpreter','latex')
xlim([xmax,0])
ylim([-0.5,2])
ytick([-0.5,0,1,2])
ylabel('$W$','Interpreter','latex')
niceplot_nobold(axfont_sz)
text(sublabelposx,1.7,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')



subplot("Position",pos(3,:));
pcolorcen(-X,Y,cFbr(:,:,tind));
xlabel('$x \, (\mathrm{m})$','Interpreter','latex')
ylabel('$y \, (\mathrm{m})$','Interpreter','latex')
xlim([xmax,0])
ylim([0,250])
title('$\nabla \times F_\mathrm{br} \,(\mathrm{m}^{-2})$','Interpreter','latex')
col=colorbar;
cmocean('balance');
%shading flat
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-0.2,0.2])
col.Ticks = linspace(-0.2,0.2,5);
xtick(-200:100:0) %modified
ytick(0:50:250)
niceplot_nobold(axfont_sz)
text(sublabelposx,230,subfiglabel{3},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
set(gca,'Layer','top')


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf -painters '/data1/bliu/vortpaper/fig_example_forcing.pdf'