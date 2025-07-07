% Bingchen Liu Feb 25, 2025
% This code analyze the CXT data with dx width and for further plotting
clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/runnum_72run.mat')
load('/data1/bliu/data/cxt_dxwidth') %'get_cxt_dxwidth_analysis'
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced.mat')

%% 5 loc with fit 
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

dx_width = -4:1:4;
trans = 0.2;

figure()
subplot(121)
for runind=1:24
    for xind = 1:5
    plot(dx_width,cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    end  

end 
hold off
xlim([0,4])
xlabel('dx lag')
ylabel('C_{X}')
niceplot_nobold(15)
colormap(cmocean('thermal',5))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.Label.FontSize = 12;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Label.Rotation = 0;



subplot(122)
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans])
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans])
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',[col(xind,:),trans])
    end  
end 
errorbar(binmean2.dx_bincenter,binmean2.cxt_binmean,binmean2.error_binmean,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m','color','m','LineStyle','-.')
hold on 
plot(fit.dx_fit,fit.cxt_fit,'LineWidth',2,'LineStyle','--','Color','k')
hold off
xlim([0,2.2])
xlabel('$dx/(c \tau)$','Interpreter','latex')
ylabel('C_{X}')
title('cxt_{fit}=exp(-x/a)*cos(x/b+c)/cos(c)')
niceplot_nobold(15)

%% 1m res (too many pts)
col = cmocean('thermal',100);
col_ind_xnond = linspace(-1,0,100);

figure()
subplot(1,3,1) % dimensional plot 
for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp2(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp2(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot([-2:1:2], cxt_dxwidth.slp2(runnum).cxt{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp3(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp3(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot([-2:1:2], cxt_dxwidth.slp3(runnum).cxt{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp4(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp4(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot([-2:1:2], cxt_dxwidth.slp4(runnum).cxt{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

hold off
clear colormap
col=colorbar;
cmocean('thermal')
caxis([-1,0])
col.Label.Interpreter = 'latex';
col.Label.String = '$ x/x_{b}$';
col.Label.FontSize = 20;
col.Label.FontWeight = 'bold';
xlabel('dx lag')
ylabel('CXT')
niceplot_nobold(15)


%
subplot(1,3,2) %nond plot 
col = cmocean('thermal',100);
col_ind_xnond = linspace(-1,0,100);

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp2(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp2(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp2(runnum).dx_nond1{1,xloc},cxt_dxwidth.slp2(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp3(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp3(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp3(runnum).dx_nond1{1,xloc},cxt_dxwidth.slp3(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp4(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp4(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp4(runnum).dx_nond1{1,xloc},cxt_dxwidth.slp4(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

hold off
clear colormap
col=colorbar;
cmocean('thermal')
caxis([-1,0])
col.Label.Interpreter = 'latex';
col.Label.String = '$ x/x_{b}$';
col.Label.FontSize = 20;
col.Label.FontWeight = 'bold';
xlabel('$dx/h$','Interpreter','latex')
ylabel('$C_{XT}/C_{XT}(dx=0)$','Interpreter','latex')
ylim([-inf,inf])
xlim([0,10])
%xlim([0,1.2])
niceplot_nobold(15)

%
subplot(1,3,3) %nond plot 
col = cmocean('thermal',100);
col_ind_xnond = linspace(-1,0,100);

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp2(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp2(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp2(runnum).dx_nond2{1,xloc},cxt_dxwidth.slp2(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp3(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp3(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp3(runnum).dx_nond2{1,xloc},cxt_dxwidth.slp3(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

for runnum=1:24
    for xloc = 1:10:length(cxt_dxwidth.slp4(runnum).dx_nond1)
        for tlag = 1:1
            x_nond = fitpara.slp4(runnum).x_nond(xloc);
            [~,col_ind] = min(abs(x_nond-col_ind_xnond));
            plot(cxt_dxwidth.slp4(runnum).dx_nond2{1,xloc},cxt_dxwidth.slp4(runnum).cxt_nond{1,xloc},'Color',col(col_ind,:))
            hold on
        end 
    end 
end 

hold off
clear colormap
col=colorbar;
cmocean('thermal')
caxis([-1,0])
col.Label.Interpreter = 'latex';
col.Label.String = '$ x/x_{b}$';
col.Label.FontSize = 20;
col.Label.FontWeight = 'bold';
xlabel('$dx/c t_{cst}$','Interpreter','latex')
ylabel('$C_{XT}/C_{XT}(dx=0)$','Interpreter','latex')
ylim([-inf,inf])
xlim([0,1])

%xlim([0,1.2])
niceplot_nobold(15)


%% 5 locations 

x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

dx_width = -4:1:4;

figure()
subplot(131)
for runind=1:24
    for xind = 1:5
    plot(dx_width,cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    end  

end 
hold off
xlim([0,4])
xlabel('dx lag')
ylabel('CXT')
niceplot_nobold(15)
colormap(cmocean('thermal',5))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.Label.FontSize = 12;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';
cbar.Label.Rotation = 0;

subplot(132)
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    end  
end 
errorbar(binmean1.dx_bincenter,binmean1.cxt_binmean,binmean1.error_binmean,'s','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m','color','m')
hold on 
errorbar(0,1,0,'s','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m','color','m')

hold off
xlim([0,1])
xlabel('dx/(c t_{cst})')
ylabel('CXT')
niceplot_nobold(15)

subplot(133)
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,1}(xind,:),'Color',col(xind,:))
    end  
end 
errorbar(binmean2.dx_bincenter,binmean2.cxt_binmean,binmean2.error_binmean,'s','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m','color','m')
hold on 
errorbar(0,1,0,'s','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m','color','m')

hold off
xlim([0,2.2])
xlabel('$dx/(c \tau)$','Interpreter','latex')
ylabel('CXT')
niceplot_nobold(15)

%% different dt lag (not working)


tlag_length  = 3;
col = cmocean('thermal',tlag_length);

cxt_dxwidt_mean = mean(cxt_dxwidth_all,1,'omitmissing');
cxt_dxwidt_std = std(cxt_dxwidth_all,1,'omitmissing');

dx_width = -4:1:4;

figure()
subplot(131)
for tlag = 1:3
for runind=1:24
    for xind = 1:5
    plot(dx_width,cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(dx_width,cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    end  
    %errorbar(dx_width(3:5),cxt_dxwidt_mean,cxt_dxwidt_std,'s','LineWidth',2.5,'MarkerSize',5,'MarkerFaceColor','m','color','m')


end 
end 
hold off
xlim([0,4])
xlabel('dx lag')
ylabel('CXT')
niceplot_nobold(15)
colormap(cmocean('thermal',3))
cbar=colorbar('south');
clim([0,2]);
cbar.Ticks =0:1:2;
cbar.Label.FontSize = 16;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$dt$';
cbar.Label.Rotation = 0;

subplot(132)
for tlag = 1:3
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond1(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    end  

end 
end 
hold off
xlim([0,1.2])
xlabel('dx/(c t_{cst})')
ylabel('CXT')
niceplot_nobold(15)

subplot(133)
for tlag=1:3
for runind=1:24
    for xind = 1:5
    plot(cxt_dxwidth_5loc.slp2(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp2(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp3(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp3(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    hold on 
    plot(cxt_dxwidth_5loc.slp4(runind).dx_nond2(xind,:),cxt_dxwidth_5loc.slp4(runind).cxt_nond{1,tlag}(xind,:),'Color',col(tlag,:))
    end  
end 
end 
hold off
xlim([0,2.2])
xlabel('$dx/(c \tau)$','Interpreter','latex')
ylabel('CXT')
niceplot_nobold(15)