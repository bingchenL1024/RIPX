% Bingchen Liu May 12, 2025
% This code compare nond tau (decorrelation time scale)

clear
close all
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'
runnum_diffslp = load('runnum_slp.mat');

[tau_nond_stats.x_mean,tau_nond_stats.binmean,tau_nond_stats.binstd] = bin_mean_taunondVSxnond(x_nond_tot_5loc,t_decoscal_tot_5loc,t_scale2_tot_5loc);

%%

xsize=13;ysize=6.8;

x0=0.1; 
y0=0.17;  
dx_sub=0.039;
xw=0.42; 
yw=0.75;
x1= x0+xw+dx_sub; 
pos = [x0 y0 xw yw; x1 y0 xw yw];

fig_fontsz= 10;
subfiglabel_fontsz = 13;
cbar_fontsize =10;
cbar_ticksz= 8;
font_rsq = 10;
subfiglabel = {'(a)','(b)'};

markersz= 15;
%%
x_num = 5;

figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

sub1= subplot("Position",pos(1,:));
for runind = 24:-1:1
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).x_nond(xind),fitpara_5loc.slp2(runind).a(xind)/fitpara_5loc.slp2(runind).t_scale2(xind),markersz,fitpara_5loc.slp2(runind).dirspr_b(xind),'filled','Marker','o')
        hold on 
        scatter(fitpara_5loc.slp3(runind).x_nond(xind),fitpara_5loc.slp3(runind).a(xind)/fitpara_5loc.slp3(runind).t_scale2(xind),markersz,fitpara_5loc.slp3(runind).dirspr_b(xind),'filled','Marker','square')
        hold on 
        scatter(fitpara_5loc.slp4(runind).x_nond(xind),fitpara_5loc.slp4(runind).a(xind)/fitpara_5loc.slp4(runind).t_scale2(xind),markersz,fitpara_5loc.slp4(runind).dirspr_b(xind),'filled','Marker','^')
        hold on 
    end 
end 
errorbar(tau_nond_stats.x_mean,tau_nond_stats.binmean,tau_nond_stats.binstd,'.','LineWidth',1.2,'Color','m')
%xlim([0,1])
%ylim([0,8])
yticks(0:2:8)
hold off
xlabel('$x/L_{sz}$','Interpreter','latex')
ylabel('$\tilde{\tau} = \hat{\tau}  / \sqrt{H_{\mathrm{sb}}/g}$','Interpreter','latex')
clear col
cbar = colormap(cmocean('thermal'));
cbar=colorbar('south');
%clim([-0.75-0.125/2,-0.25+0.125/2]);
%cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.FontSize= cbar_ticksz;
cbar.Ruler.TickLabelRotation=0;
%cbar.FontAngle = 0; 
cbar.Label.FontSize = cbar_fontsize;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\sigma_{\theta b} \, (\mathrm{deg})$';
cbar.Position = [0.17, 0.87,0.28,0.035];
grid on 
%title('$\alpha = 0.5$','Interpreter','latex')
niceplot_nobold(fig_fontsz)
text(-0.78,7.5,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
clear col




col = cmocean('thermal',x_num);
sub2= subplot("Position",pos(2,:));
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).G0_nond(xind),fitpara_5loc.slp2(runind).a(xind)/fitpara_5loc.slp2(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','o')
        hold on 
        scatter(fitpara_5loc.slp3(runind).G0_nond(xind),fitpara_5loc.slp3(runind).a(xind)/fitpara_5loc.slp3(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','square')
        hold on 
        scatter(fitpara_5loc.slp4(runind).G0_nond(xind),fitpara_5loc.slp4(runind).a(xind)/fitpara_5loc.slp4(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','^')
        hold on 
    end 
end 
plot(tau_nond_fit.G0_nond_fit,tau_nond_fit.tau_nond_fit,'LineWidth',3,'Color','r','LineStyle','--')
xlim([0,1])
%ylim([0,8])
%xticks(linspace(1,7,7))
%yticks(linspace(1,7,7))
hold off
yticks(0:2:8)
xlabel('$\hat{G}_{0}/\hat{G}_{\mathrm{max}}$','Interpreter','latex')
%ylabel('$\tilde{\tau} = \hat{\tau} \sqrt{g/h}$','Interpreter','latex')
yticklabels([])
clear col
cbar2=colormap(sub2,cmocean('thermal',x_num));
cbar2=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar2.Ticks =linspace(-0.75,-0.25,5);
cbar2.FontSize= cbar_ticksz;
cbar2.Ruler.TickLabelRotation=0;
%cbar2.FontAngle = 0; 
cbar2.Label.FontSize = cbar_fontsize;
cbar2.Label.FontWeight = 'bold';
cbar2.Label.Interpreter = 'latex';
cbar2.Label.String = '$x / L_{\mathrm{sz}}$';
cbar2.Position = [0.63, 0.87,0.28,0.035];
grid on 
niceplot_nobold(fig_fontsz)
text(0.035,7.5,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')



set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf -painters '/data1/bliu/vortpaper/fig_tau_nond.pdf'


