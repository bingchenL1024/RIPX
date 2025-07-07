% Bingchen Liu March 11, 2025
% This code compare the linear regression of dx location of max CXT to
% c=\sqrt(gh) prediction to validate the use of choice of max CXT 

clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth') %'get_cxt_alongct_nointerp_max_dxwidth'
load('/data1/bliu/data/runnum_72run')

RMSE_dist = Stats.RMSE_dist
RMSE_vel = Stats.RMSE_vel
rsquare = Stats.c_corr
%%

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

markersz= 20;
%%
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(1,:));
for xloc =x_num:-1:1
    for N =1:24
        scatter(c_phase_5loc.slp2(N).c_modelh(xloc),c_phase_5loc.slp2(N).c_fit(xloc),markersz,col(xloc,:),'filled','MarkerEdgeColor','k','Marker','o')
        hold on 
        scatter(c_phase_5loc.slp3(N).c_modelh(xloc),c_phase_5loc.slp3(N).c_fit(xloc),markersz,col(xloc,:),'filled','MarkerEdgeColor','k','Marker','square')
        hold on 
        scatter(c_phase_5loc.slp4(N).c_modelh(xloc),c_phase_5loc.slp4(N).c_fit(xloc),markersz,col(xloc,:),'filled','MarkerEdgeColor','k','Marker','^')
        hold on 
    end
end 
xlim([0,7])
ylim([0,7])
xticks(linspace(1,7,7))
yticks(linspace(1,7,7))
plot_oneone
hold off
xlabel('$\sqrt{g h} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
ylabel('$\hat{c} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
clear col
colormap(cmocean('thermal',x_num))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.FontSize= cbar_ticksz;
cbar.Ruler.TickLabelRotation=0;
%cbar.FontAngle = 0; 
cbar.Label.FontSize = cbar_fontsize;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Position = [0.38, 0.18,0.39,0.02];
grid on 
niceplot_nobold(fig_fontsz)


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf -painters '/data1/bliu/vortpaper/fig_cxt_c_compare.pdf'






