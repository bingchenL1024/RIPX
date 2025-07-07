% Bingchen Liu May 12, 2025
% This code compare nond tau (decorrelation time scale)

clear
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'


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
font_rsq = 10;

markersz= 15;
%%
x_num = 5;
col = cmocean('thermal',x_num);

figure()

set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(1,:));
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).G0_nond(xind),fitpara_5loc.slp2(runind).a(xind)/fitpara_5loc.slp2(runind).t_scale(xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp3(runind).G0_nond(xind),fitpara_5loc.slp3(runind).a(xind)/fitpara_5loc.slp3(runind).t_scale(xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp4(runind).G0_nond(xind),fitpara_5loc.slp4(runind).a(xind)/fitpara_5loc.slp4(runind).t_scale(xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
        hold on 
    end 
end 
plot(tau_nond_fit.G0_nond_fit,tau_nond_fit.tau_nond_fit,'LineWidth',3,'Color','r','LineStyle','--')
xlim([0,1])
%ylim([0,7])
%xticks(linspace(1,7,7))
%yticks(linspace(1,7,7))

hold off
xlabel('$G_{0}/G_{\mathrm{max}}$','Interpreter','latex')
ylabel('$\tilde{\tau} = \hat{\tau} \sqrt{g/h}$','Interpreter','latex')
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
cbar.Position = [0.38, 0.85,0.39,0.02];
grid on 
niceplot_nobold(fig_fontsz)


set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf -painters '/data1/bliu/vortpaper/fig_tau_nond.pdf'


