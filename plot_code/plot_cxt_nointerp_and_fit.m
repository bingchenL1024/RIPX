% Bingchen Liu Nov 13, 2024 
% This code plot cxt without interpolation and its fit for visualization 

clear

%% dimensional Cxt without fit with dt 
load('/data1/bliu/data/cxt_alongct_nointerp')
load('/data1/bliu/data/cxt_ind_good.mat')


runnum = 1;
skipx = 3;

ind_good = ind_good_All{runnum};
cxt_alongct = cxt_alongct_ALL{runnum};
x= x_nond_All{runnum};
x_cxt = x_cxt_alongct_All{runnum};

colormap = cmocean('solar',length(ind_good)); %parula

figure()
subplot(121)
for i = 1:skipx:length(ind_good)
    plot(t_itp, cxt_alongct(:,i),'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
clear colormap
colormap(cmocean('solar',length(ind_good)));
%legend([p1, p2, p3],'Mean', '25 percentile','75 percentile')
xlabel('$\Delta t (s)$','Interpreter','latex')
ylabel('$\tilde{C}_{XT}$','Interpreter','latex')
%title(['Run',num2str(runnum)])
colorbar
col=colorbar;
col.Label.FontSize = 22;
caxis([x(1),x(end)])
col.Label.Interpreter = 'latex';
col.Label.String = 'Non-dimensional Cross-shore Location $\frac{x}{x_{b}}$';
%xlim([0,9])
niceplot_nobold_nomintick(18);
grid on 
hold off 

colormap = cmocean('solar',length(ind_good)); %parula
subplot(122)
for i = 1:skipx:length(ind_good)
    plot(x_cxt(:,i), cxt_alongct(:,i),'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
clear colormap
colormap(cmocean('solar',length(ind_good)));
%legend([p1, p2, p3],'Mean', '25 percentile','75 percentile')
xlabel('$\Delta x (m)$','Interpreter','latex')
ylabel('$\tilde{C}_{XT}$','Interpreter','latex')
%title(['Run',num2str(runnum)])
colorbar
col=colorbar;
col.Label.FontSize = 22;
caxis([x(1),x(end)])
col.Label.Interpreter = 'latex';
col.Label.String = 'Non-dimensional Cross-shore Location $\frac{x}{x_{b}}$';
%xlim([0,9])
niceplot_nobold_nomintick(18);
grid on 
hold off 





%% Cxt WITH fit 


load('/data1/bliu/data/cxt_alongct_nointerp_fitpara_qced.mat')
load('/data1/bliu/data/cxt_ind_good.mat')

runnum = 5;
xind_skip = 3;

% ind_good = ind_good_All{runnum};
% cxt_alongct = cxt_alongct_ALL{runnum};
% x= x_nond_All{runnum};
% x_cxt = x_cxt_alongct_All{runnum};

colormap = cmocean('solar',length(fitpara.slp2(runnum).cxt_data)); %parula

figure()
for xloc = 1:xind_skip:length(fitpara.slp2(runnum).cxt_data)
    plot(fitpara.slp2(runnum).t_itp{xloc}, fitpara.slp2(runnum).cxt_data{xloc},'LineWidth',1.5,'Color',colormap(xloc,:))
    hold on 
    plot(fitpara.slp2(runnum).t_itp{xloc},fitpara.slp2(runnum).cxt_fit{xloc},'--','LineWidth',1.5,'Color',colormap(xloc,:))
    hold on 
    legend('Data','AR1 Model Fit')
end 
clear colormap
colormap(cmocean('solar',length(ind_good)));
%legend([p1, p2, p3],'Mean', '25 percentile','75 percentile')
xlabel('$\Delta t (s)$','Interpreter','latex')
ylabel('$\tilde{C}_{XT}$','Interpreter','latex')
%title(['Run',num2str(runnum)])
colorbar
col=colorbar;
col.Label.FontSize = 22;
caxis([fitpara.slp2(runnum).x_nond(1),fitpara.slp2(runnum).x_nond(end)])
col.Label.Interpreter = 'latex';
col.Label.String = 'Non-dimensional Cross-shore Location $\frac{x}{x_{b}}$';
%xlim([0,9])
niceplot_nobold_nomintick(18);
grid on 
hold off 