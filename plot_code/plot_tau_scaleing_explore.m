%Bingchen Liu july 4, 2025
%this code test the dimensional scaling for \tau

clear
close all
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'

%% tau VS scaling

xsize=8;ysize=8;
x0=0.15; 
y0=0.15;  
dx=0.12;
xw=0.75; yw=xw*(xsize/ysize);
x1= x0+xw+dx; 
pos = [x0 y0 xw yw];

fig_fontsz= 10;
subfiglabel_fontsz = 15;
cbar_fontsize =10;
cbar_ticksz= 6;
font_rsq = 10;

markersz= 35;

x_num = 5;
col = cmocean('thermal',x_num);

figure()
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).t_scale2(xind),fitpara_5loc.slp2(runind).a(xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp3(runind).t_scale2(xind),fitpara_5loc.slp3(runind).a(xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp4(runind).t_scale2(xind),fitpara_5loc.slp4(runind).a(xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
        hold on 
    end 
end 
ylabel('$\tau$','Interpreter','latex')
%xlabel('$g h^2 D_{w}^{-1}$','Interpreter','latex')
xlabel('$D_{w}^{1/3}/g$','Interpreter','latex')
%xlim([0,1])
%ylim([0,7])
%xticks(linspace(1,7,7))
%yticks(linspace(1,7,7))

hold off