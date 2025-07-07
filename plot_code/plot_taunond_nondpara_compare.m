% Bingchen Liu July 4, 2025
% This code compare nond tau (decorrelation time scale) VS nond parameters
% of
% Dw/Dw_max
% cFbr/cFbr_max where cFbr is parameterized using Dw and paper scaling 

clear
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc') %'get_cxt_5loc'
%% test figure(quick, no detailed color and symbol)
% figure()
% scatter(t_scale2_tot_5loc,t_decoscal_tot_5loc)
% 
% 


%% tau_nond VS nond parameter

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

%subplot(211)
for runind = 1:24
    for xind = 1:5
        scatter(fitpara_5loc.slp2(runind).x_nond(xind),fitpara_5loc.slp2(runind).a(xind).*0.02/fitpara_5loc.slp2(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp3(runind).x_nond(xind),fitpara_5loc.slp3(runind).a(xind).*0.03/fitpara_5loc.slp3(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
        hold on 
        scatter(fitpara_5loc.slp4(runind).x_nond(xind),fitpara_5loc.slp4(runind).a(xind).*0.04/fitpara_5loc.slp4(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
        hold on 
    end 
end 
%xlim([0,1])
%xticks(linspace(1,7,7))
%yticks(linspace(1,7,7))

hold off
xlabel('$x/x_b$','Interpreter','latex')
%ylabel('$\tilde{\tau} = \hat{\tau} / gh^2 /D_w$','Interpreter','latex')
%ylabel('$\tilde{\tau} = \hat{\tau} \beta / D_{w}^{1/3} g$','Interpreter','latex')
%ylabel('$\tilde{\tau} = \hat{\tau} / \sigma_{\theta}/\sigma_{b}$','Interpreter','latex')

clear col
colormap(cmocean('thermal',x_num))
cbar=colorbar;
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.FontSize= cbar_ticksz;
cbar.Ruler.TickLabelRotation=0;
%cbar.FontAngle = 0; 
cbar.Label.FontSize = cbar_fontsize;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
%cbar.Position = [0.38, 0.85,0.39,0.02];
grid on 
niceplot_nobold(fig_fontsz)

% col = cmocean('thermal',x_num);
% subplot(212)
% for runind = 1:24
%     for xind = 1:5
%         scatter(fitpara_5loc.slp2(runind).curlFbr_para_nond(xind),fitpara_5loc.slp2(runind).a(xind)/fitpara_5loc.slp2(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
%         hold on 
%         scatter(fitpara_5loc.slp3(runind).curlFbr_para_nond(xind),fitpara_5loc.slp3(runind).a(xind)/fitpara_5loc.slp3(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
%         hold on 
%         scatter(fitpara_5loc.slp4(runind).curlFbr_para_nond(xind),fitpara_5loc.slp4(runind).a(xind)/fitpara_5loc.slp4(runind).t_scale2(xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
%         hold on 
%     end 
% end 
% xlim([0,1])
% hold off
% xlabel('$\mathrm{cFbr}/\mathrm{cFbr}_\mathrm{max}$','Interpreter','latex')
% %ylabel('$\tilde{\tau} = \hat{\tau} / gh^2 /D_w$','Interpreter','latex')
% ylabel('$\tilde{\tau} = \hat{\tau} / h D_{w}^{-1/3}$','Interpreter','latex')
% clear col
% colormap(cmocean('thermal',x_num))
% cbar=colorbar;
% clim([-0.75-0.125/2,-0.25+0.125/2]);
% cbar.Ticks =linspace(-0.75,-0.25,5);
% cbar.FontSize= cbar_ticksz;
% cbar.Ruler.TickLabelRotation=0;
% %cbar.FontAngle = 0; 
% cbar.Label.FontSize = cbar_fontsize;
% cbar.Label.FontWeight = 'bold';
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$x / L_{\mathrm{sz}}$';
% %cbar.Position = [0.38, 0.85,0.39,0.02];
% grid on 
% title('$\alpha = 0.5$','Interpreter','latex')
% niceplot_nobold(fig_fontsz)