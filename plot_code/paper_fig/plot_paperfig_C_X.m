% Bingchen Liu Mar 27, 2025
% this code plot C_X part for paper 

clear
close all 
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/runnum_72run.mat')
load('/data1/bliu/data/cxt_dxwidth') %'get_cxt_dxwidth_analysis'
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced.mat')

%% 5 loc with fit 
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

dx_width = -4:1:4;
trans = 0.2;
xsize=13;ysize=6.8;

subfiglabel = {'(a)','(b)'};
linewidth = 0.4;
linewidth_thick =1;
fig_fontsize = 10;
subfiglabel_fontsz = 12;
leg_fontsz = 8;

x0=0.13; 
y0=0.17;  
dx_sub=0.03;
xw=0.41; 
yw=0.75;
x1= x0+xw+dx_sub; 
pos = [x0 y0 xw yw; x1 y0 xw yw];

%============================  plot ================
close all


figure()
set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

subplot("Position",pos(1,:));
for runind=1:24
    for xind = 1:5
    plot(dx_width,cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:),'LineWidth',linewidth)
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:),'LineWidth',linewidth)
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:),'LineWidth',linewidth)
    end  

end 
hold off
grid on 
xlim([0,4])
ylim([-0.8,1])
set(gca,'YTick',[-0.5:0.5:1])
set(gca,'XTick',[0:1:4])
xlabel('$\Delta x''\, (\mathrm{m})$','Interpreter','latex')
ylabel('$C_{\textit{X}}$','Interpreter','latex')
niceplot_nobold(fig_fontsize)
text(3.5,0.91,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')
colormap(cmocean('thermal',5))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
%cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Label.Rotation = 0;
cbar.Position = [0.18, 0.88,0.29,0.03];
cbar.Label.Rotation = 0;
cbar.Ruler.TickLabelRotation=0;
cbar.FontSize = 7;
cbar.Label.FontSize = 12;



subplot("Position",pos(2,:));
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans],'LineWidth',linewidth)
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans],'LineWidth',linewidth)
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans],'LineWidth',linewidth)
    end  
end 
errorbar(binmean2.dx_bincenter,binmean2.cxt_binmean,binmean2.error_binmean,'-s','LineWidth',linewidth_thick-0.3,'MarkerSize',3,'MarkerFaceColor','m','color','m','LineStyle','-.')
hold on 
plot(fit.dx_fit,fit.cxt_fit,'LineWidth',linewidth_thick,'LineStyle','--','Color','r')
hold off
grid on 
xlim([0,2.2])
ylim([-0.8,1])
set(gca,'YTick',[-0.5:0.5:1])
set(gca,'XTick',[0:0.5:2])
xlabel('$\Delta x''/(\hat{c} \, \hat{\tau})$','Interpreter','latex')
ylabel('')
yticklabels([])
niceplot_nobold(fig_fontsize)
text(1.9,0.91,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')



set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_C_X_5loc.pdf'

close all














