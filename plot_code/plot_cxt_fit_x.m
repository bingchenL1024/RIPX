% Bingchen Liu, Aug 14, 2024
% This code plot the cxt and its fit parameters
% and evaluate the fit parameters with dimensionless number 
clear
load('/data1/bliu/data/cxt_alongct_x_fitpara_qced.mat')

%%
figure()
for run = 1:4:24
plot(fitpara.slp2(run).x,fitpara.slp2(run).a)
hold on 
end 
hold off 
%% cxt(dimensional x) amd include cxt_fit 

runnum = 5;
xind_skip = 5;

clear colormap
colormap = colormap(cmocean('haline',length(fitpara.slp2(runnum).x_itp_fit)));
close all
figure()
for xloc = 1:xind_skip:length(fitpara.slp2(runnum).x_itp_fit)
    plot(fitpara.slp2(runnum).x_itp_fit{xloc},fitpara.slp2(runnum).cxt_data{xloc},'LineWidth',1.5,'Color',colormap(xloc,:))
    hold on 
    plot(fitpara.slp2(runnum).x_itp_fit{xloc},fitpara.slp2(runnum).cxt_fit{xloc},'--','LineWidth',1.5,'Color',colormap(xloc,:))
    hold on 
    legend('Data','AR2 Model Fit')
end 
clear colormap
colormap(cmocean('halin',length(fitpara.slp2(runnum).x_itp_fit)));
xlabel('$\Delta x (m)$','Interpreter','latex')
ylabel('$\tilde{C}_{XT}$','Interpreter','latex')
%title('$cxt_{fit} = \exp{(-\frac{x}{a})} \cos{(\frac{x}{b} +c)}/\cos{c}$','interpreter','latex')
colorbar
col=colorbar;
col.Label.FontSize = 22;
caxis([fitpara.slp2(runnum).x_nond(1),fitpara.slp2(runnum).x_nond(end)])
col.Label.Interpreter = 'latex';
col.Label.String = 'Non-dimensional Cross-shore Location $\frac{x}{x_{b}}$';niceplot_nobold_nomintick(18);
grid on 
hold off 
niceplot(18)

%%  cxt(dimensional x)

load('/data1/bliu/data/cxt_alongct_x.mat')

runnum = 1;
skipx = 1;

ind_good = ind_good_All{runnum};
cxt_alongct = cxt_alongct_x_ALL{runnum};
x= x_nond_All{runnum};

colormap = cmocean('solar',length(ind_good)); %parula

figure()
for i = 1:skipx:length(ind_good)
    plot(x_itp, cxt_alongct(:,i),'LineWidth',1.5,'Color',colormap(i,:))
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


%% cxt(nond x)

load('/data1/bliu/data/cxt_alongct_x.mat')
load('/data1/bliu/data/cxt_alongct_x_fitpara','fitpara','readme')

runnum =1;
nond_x = fitpara(runnum).slp2.x_nond;
x_itp = 0:0.1:9;
kw = fitpara(runnum).slp2.kw;
h = fitpara(runnum).slp2.h;

cxt = fitpara(runnum).slp2.cxt_data;
for i  = 1:length(kw)
    nond_dx(:,i) = x_itp.*kw(i);
    %nond_dx(:,i) = x_itp./h(i);
end 

colormap = cmocean('solar',length(nond_x)); %parula

figure()
for i = 1:length(nond_x)
    plot(nond_dx(:,i), cxt{i},'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
clear colormap
colormap(cmocean('solar',length(nond_x)));
%legend([p1, p2, p3],'Mean', '25 percentile','75 percentile')
xlabel('$\Delta x \, k_{w}$','Interpreter','latex')
ylabel('$\tilde{C}_{XT}$','Interpreter','latex')
%title(['Run',num2str(runnum)])
colorbar
col=colorbar;
col.Label.FontSize = 22;
caxis([nond_x(1),nond_x(end)])
col.Label.Interpreter = 'latex';
col.Label.String = 'Non-dimensional Cross-shore Location $\frac{x}{x_{b}}$';
%xlim([0,9])
niceplot_nobold_nomintick(18);
grid on 
hold off 
