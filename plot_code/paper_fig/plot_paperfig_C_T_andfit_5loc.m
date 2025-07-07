% Bingchen Liu Mar 12, 2025
% This code plot the cxt function and its fit for paper


clear
load('/data1/bliu/data/cxt_alongct_nointerp')
load('/data1/bliu/data/cxt_ind_good.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc')
load('/data1/bliu/data/ind_of_diff_bath.mat')

%%
%co = [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4660 0.6740 0.1880]; % blue-yellow-green
co = cmocean('thermal',5);
subfiglabel = {'(a)','(b)'};
linewidth = 0.5;
fig_fontsize = 10;
subfiglabel_fontsz = 12;
leg_fontsz = 8;
trans=0.2;
t_nond_fit = 0:0.1:12;
linewidth_thick =1;


xsize=12;ysize=6.8;

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

% ================================= sub 2 =================================
%subplot(121)
subplot("Position",pos(1,:));
for i = 1:24
    for loc = 1:5
        h1=plot(fitpara_5loc.slp2(i).t_itp{loc},fitpara_5loc.slp2(i).cxt_data{loc},'Color',co(loc,:),'LineWidth',linewidth);
        hold on 
        h2=plot(fitpara_5loc.slp3(i).t_itp{loc},fitpara_5loc.slp3(i).cxt_data{loc},'Color',co(loc,:),'LineWidth',linewidth);
        hold on 
        h3=plot(fitpara_5loc.slp4(i).t_itp{loc},fitpara_5loc.slp4(i).cxt_data{loc},'Color',co(loc,:),'LineWidth',linewidth);
        hold on 
    end 
end 
xlim([0,5])
ylim([0,1])
ylabel('$C_{\textit{T}}$','Interpreter','latex')
xlabel('$\Delta t \, (s)$','Interpreter','latex')
xtickangle(35)
set(gca,'XTick',0:1:5)
set(gca,'YTick',[0:0.2:1])
ax=gca;
set(ax.YLabel,'Position',[-0.9,ax.YLabel.Position(2),0])
%set(ax.XLabel,'Position',[ax.XLabel.Position(1),-0.3,0])
xtickangle(0)
grid on
hold off 
niceplot_nobold_nomintick(fig_fontsize)
clear col
colormap(cmocean('thermal',5))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ruler.TickLabelRotation=0;
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.FontSize = 7;
cbar.Label.FontSize = 12;
%cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Label.Rotation = 0;
cbar.Position = [0.18, 0.88,0.29,0.03];
text(4.3,0.93,subfiglabel{1},'FontSize',subfiglabel_fontsz,'Interpreter','latex')

% ================================= sub 2 =================================
%subplot(122)
co = cmocean('thermal',5);

subplot("Position",pos(2,:));
for i = 1:24
    for loc =1:5
    plot(fitpara_5loc.slp2(i).t_nond_diva{loc},fitpara_5loc.slp2(i).cxt_data{loc},'Color',[co(loc,:),trans],'LineWidth',linewidth-0.1)
    hold on 

    plot(fitpara_5loc.slp3(i).t_nond_diva{loc},fitpara_5loc.slp3(i).cxt_data{loc},'Color',[co(loc,:),trans],'LineWidth',linewidth-0.1)
    hold on 

    plot(fitpara_5loc.slp4(i).t_nond_diva{loc},fitpara_5loc.slp4(i).cxt_data{loc},'Color',[co(loc,:),trans],'LineWidth',linewidth-0.1)
    hold on 
    end 

end 
plot(t_nond_fit,exp(-t_nond_fit),'LineWidth',linewidth_thick,'LineStyle','--','Color','r')

xlim([0,8])
ylim([0,1])
ylabel('')
xlabel('$\Delta t / \hat{\tau}$','Interpreter','latex')
yticklabels([])
set(gca,'XTick',0:2:12)
set(gca,'YTick',[-0.1,0:0.2:1])
ax=gca;
set(ax.YLabel,'Position',[-2.7,ax.YLabel.Position(2),0])
grid on 
hold off 
xtickangle(0)
niceplot_nobold_nomintick(fig_fontsize)
text(10.0,0.93,subfiglabel{2},'FontSize',subfiglabel_fontsz,'Interpreter','latex')






set(gcf,'PaperUnit','centimeters')
set(gcf,'papersize',[xsize,ysize],'PaperPosition',[0 0 xsize ysize])

print -dpdf '/data1/bliu/vortpaper/fig_C_T_5loc.pdf'

close all
