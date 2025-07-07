% Bingchen Liu March 12, 2025
% This code plot the 4 panels CXT(dt,dx) plot for paper fig
clear
load('/data1/bliu/data/cxt_atx_4panels') %generated in get_4panel_cxtatx


xsize=15;ysize=15;

x0=0.11;  
dx=0.03;
xw=0.405;
x1 = x0+xw+dx;
y0=0.1;  
yw=0.38; 
dy=0.08;
y1 = y0+dy+yw;  y2 = y1+dy+yw;
pos = [x0 y1 xw yw; x0 y0 xw yw; x1 y1 xw yw; x1 y0 xw yw];

label_locx = 0.15;
label_locy = 14;

subfiglabel = {'(a)','(b)','(c)','(d)'};


fontsz = 12;
subfiglabel_fontsz = 18;
xlim1=7;
ylim1=15;
x_lag = 0:20; % in meters
t_lag = -10:10; %in seconds
linewidth_c = 1.2;

[dt,dx] = meshgrid(t_lag,x_lag);



figure()
ind=1;
subplot("Position",pos(1,:))
ct = c(ind).*t_lag;
ct_fit = c_fit_4loc(ind).*t_lag;
pcolorcen(dt,dx,cxt_atx(:,:,ind));
hold on 
plot(t_lag,ct,'LineWidth',linewidth_c,'color','k','LineStyle','--')
hold on 
plot(t_lag,ct_fit,'LineWidth',linewidth_c,'color','k','LineStyle','-.')
hold on 
scatter(dt_max_loc,dx_max_loc(:,ind),150,'rp','filled')
col=colorbar('south');
col.AxisLocation='in';
cmocean('balance');
col.Position = [0.32, 0.58,0.15,0.03];
caxis([-1,1])
col.Label.Interpreter = 'latex';
col.Label.String = '$\tilde{C}_{XT}$';
col.Label.FontSize = 14;
col.Label.FontWeight = 'bold';
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
xticks(0:2:10)
ytick(0:3:20)
ylabel('$\Delta x$ (m)','Interpreter','latex')
%xlabel('$\Delta t$ (s)','Interpreter','latex')
%title(['xlocation = ',num2str(x_nond(ind_good(i))),' dimensionless sz loc'])
ylim([0,ylim1])
xlim([0,xlim1])
set(gca,'XTickLabel',[])
title(['$x/L_{\mathrm{sz}} =$',num2str(x_nond_ind(ind))],'Interpreter','latex')
niceplot_nobold(fontsz);
text(label_locx,label_locy,subfiglabel{ind},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

subplot("Position",pos(3,:))
ind=2;
ct = c(ind).*t_lag;
ct_fit = c_fit_4loc(ind).*t_lag;
pcolorcen(dt,dx,cxt_atx(:,:,ind));
hold on 
plot(t_lag,ct,'LineWidth',linewidth_c,'color','k','LineStyle','--')
hold on 
plot(t_lag,ct_fit,'LineWidth',linewidth_c,'color','k','LineStyle','-.')
hold on 
scatter(dt_max_loc,dx_max_loc(:,ind),150,'rp','filled')
col=colorbar;
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-1,1])
col.Label.Interpreter = 'latex';
col.Label.String = '$\tilde{C}_{XT}$';
col.Label.FontSize = 22;
col.Label.FontWeight = 'bold';
colorbar off
xticks(0:2:10)
ytick(0:3:20)
%ylabel('$\Delta x$ (m)','Interpreter','latex')
%xlabel('$\Delta t$ (s)','Interpreter','latex')
%title(['xlocation = ',num2str(x_nond(ind_good(i))),' dimensionless sz loc'])
ylim([0,ylim1])
xlim([0,xlim1])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
title(['$x/L_{\mathrm{sz}} =$',num2str(x_nond_ind(ind))],'Interpreter','latex')
niceplot_nobold(fontsz);
text(label_locx,label_locy,subfiglabel{ind},'FontSize',subfiglabel_fontsz,'Interpreter','latex')


subplot("Position",pos(2,:))
ind=3;
ct = c(ind).*t_lag;
ct_fit = c_fit_4loc(ind).*t_lag;
pcolorcen(dt,dx,cxt_atx(:,:,ind));
hold on 
plot(t_lag,ct,'LineWidth',linewidth_c,'color','k','LineStyle','--')
hold on 
plot(t_lag,ct_fit,'LineWidth',linewidth_c,'color','k','LineStyle','-.')
hold on 
scatter(dt_max_loc,dx_max_loc(:,ind),150,'rp','filled')
col=colorbar;
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-1,1])
col.Label.Interpreter = 'latex';
col.Label.String = '$\tilde{C}_{XT}$';
col.Label.FontSize = 22;
col.Label.FontWeight = 'bold';
xticks(0:2:10)
ytick(0:3:20)
ylabel('$\Delta x$ (m)','Interpreter','latex')
xlabel('$\Delta t$ (s)','Interpreter','latex')
colorbar off
%title(['xlocation = ',num2str(x_nond(ind_good(i))),' dimensionless sz loc'])
ylim([0,ylim1])
xlim([0,xlim1])
title(['$x/L_{\mathrm{sz}} =$',num2str(x_nond_ind(ind))],'Interpreter','latex')
niceplot_nobold(fontsz);
text(label_locx,label_locy,subfiglabel{ind},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

subplot("Position",pos(4,:))
ind=4;
ct = c(ind).*t_lag;
ct_fit = c_fit_4loc(ind).*t_lag;
pcolorcen(dt,dx,cxt_atx(:,:,ind));
hold on 
plot(t_lag,ct,'LineWidth',linewidth_c,'color','k','LineStyle','--')
hold on 
plot(t_lag,ct_fit,'LineWidth',linewidth_c,'color','k','LineStyle','-.')
hold on 
scatter(dt_max_loc,dx_max_loc(:,ind),150,'rp','filled')
col=colorbar('south');
cmocean('balance');
%caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
caxis([-1,1])
colorbar off
xticks(0:2:10)
ytick(0:3:20)
%ylabel('$\Delta x$ (m)','Interpreter','latex')
xlabel('$\Delta t$ (s)','Interpreter','latex')
%title(['xlocation = ',num2str(x_nond(ind_good(i))),' dimensionless sz loc'])
ylim([0,ylim1])
xlim([0,xlim1])
set(gca,'YTickLabel',[])
title(['$x/L_{\mathrm{sz}} =$',num2str(x_nond_ind(ind))],'Interpreter','latex')
niceplot_nobold(fontsz);
text(label_locx,label_locy,subfiglabel{ind},'FontSize',subfiglabel_fontsz,'Interpreter','latex')



set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])


print -dpdf '/data1/bliu/vortpaper/fig_cxt_dtdx_4panels.pdf'

close all