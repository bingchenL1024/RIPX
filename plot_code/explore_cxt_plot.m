% Bingchen Liu July 15, 2024
% This is a exploratory code to plot c_xt 
clear
close all

%load ('/data1/bliu/data/raw/CXT_ALL_Mark.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_Falk.mat')
%%
load('/data1/bliu/data/SS_raw.mat')
% load data
runnum = 2;
g = 9.81;
dim = size(cell2mat(CXT_ALL(runnum)));
SS = S(runnum);
x = SS.X;
h = SS.h;
cxt = cell2mat(CXT_ALL(runnum));


cxt_min = zeros(dim(1),1);
cxt_max = zeros(dim(1),1);
for i = 1:dim(1)
    slice = cxt(i,:,:);
    cxt_max(i) = max(slice(:));
    cxt_min(i) = min(slice(:));
end 
ind_good_temp = find(cxt_max<=1&cxt_min>=-1);
ind_good_diff = diff(ind_good_temp);
ind_good = ind_good_temp(find(ind_good_diff==1));


x_lag = 0:9; % in meters
t_lag = -50:50; %in seconds
    
[dt,dx] = meshgrid(t_lag,x_lag);
%% quality check 
figure(2)
histogram(cxt(:));
xlim([-1.5,1.5]);
ylim([0,100])
xlabel('c_{xt}')
ylabel('Count')
title('Hist of cxt of all cross-shore loc (1 run)','FontSize',20)
niceplot_nobold_nomintick(22)


figure()
plot(cxt_max)
hold on 
plot(cxt_min)
ylim([-1.5,1.5])
xlabel('run index for cross-shore locations','FontSize',18)
title('Max and min of cxt at different cross shore location','FontSize',22)
hold off
%% C_XT(dx,dt) plot
for xloc=1:10:dim(1) %cross shore dim
    %xloc=  111;
    h_xloc = h(xloc);
    c = sqrt(g*h_xloc);

    cxt_atx = squeeze(cxt(xloc,:,:));
    
    ct = c.*t_lag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%cxt(dx,dt)
    contourf(dt,dx,cxt_atx,20,'LineColor','none')
    hold on 
    plot(t_lag,ct,'LineWidth',2,'color','k')
    
    
    col=colorbar;
    cmocean('balance');
    %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
    caxis([-1,1])
    col.Label.String = 'Cross Correlation';
    ylabel('x lag (m)')
    xlabel('Time lag (s)')
    title(['xlocation = ',num2str(x(xloc)),' m'])
    ylim([0,9])
    xlim([-8,8])
    niceplot_nobold_nomintick(18);
    drawnow
    pause(1)
    
end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% use pcolor
for xloc=1:10:dim(1) %cross shore dim
    %xloc=  111;
    h_xloc = h(xloc);
    c = sqrt(g*h_xloc);

    cxt_atx = squeeze(cxt(xloc,:,:));
    x_lag = 0:9; % in meters
    t_lag = -50:50; %in seconds

%     x_lag_pcolor = 0.5:9.5;
%     t_lag_pcolor = -50.5:49.5;
    
    %[dt_pcolor,dx_pcolor] = meshgrid(t_lag_pcolor,x_lag_pcolor);
    [dt_pcolor,dx_pcolor] = meshgrid(t_lag,x_lag);

    
    ct = c.*t_lag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%cxt(dx,dt)
    pcolor(dt_pcolor,dx_pcolor, cxt_atx)
    hold on 
    plot(t_lag,ct,'LineWidth',2,'color','k')
    shading interp
    
    col=colorbar;
    cmocean('balance');
    caxis([-1,1])

    
    col.Label.String = 'Cross Correlation';
    ylabel('$\Delta x$ (m)','Interpreter','latex')
    xlabel('$\Delta t$ (s)','Interpreter','latex')
    %title(['xlocation = ',num2str(x(xloc)),' m'])
    ylim([0,9])
    xlim([-8,8])
    niceplot(18);
    drawnow
    pause(1)
    
end

%% R(d\theta)
% c_xt = cell2mat(CXT_ALL(runnum));
% x_lag = 0:9; % in meters
% t_lag = -50:50; %in seconds
% 
% [dx,dt] = meshgrid(t_lag,x_lag);
% dx_vec = dx(:);
% dt_vec = dt(:);
% 
% for xloc=1:dim(1) %cross shore dim
%     h_xloc = h(xloc);
%     c(xloc) = -sqrt(g*h_xloc);
% 
%     cxt_atx = squeeze(c_xt(xloc,:,:));
%     
%     
%     cxt_atx_vec(:,xloc) = cxt_atx(:);
%     d_theta(:,xloc) = (dx_vec-c(xloc).*dt_vec)/c(xloc);
% end 
% 
% colormap = parula(dim(1));
% figure()
% for xloc=1:dim(1) %cross shore dim
% scatter(d_theta(:,xloc),cxt_atx_vec(:,xloc),'filled')
% hold on 
% colororder(colormap)
% end 
% xlim([-15,15])
% ylabel('$c_{XT}$','interpreter','latex')
% xlabel('$\Delta x - c \Delta t$','interpreter','latex')
% cb=colorbar;
% ylabel(cb,'Cross-shore location (m)')
% caxis([1,dim(1)])
% title(['run',num2str(runnum)])
% niceplot_nobold_nomintick(22);
% hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%cxt(dx,dt)
%% plot R(d/theta)
% for xloc=1:10:dim(1) %cross shore dim
% scatter(d_theta(:,xloc),cxt_atx_vec(:,xloc),'filled')
% 
% xlim([-30,0])
% ylabel('$c_{XT}$','interpreter','latex')
% xlabel('$\Delta x - c \Delta t$','interpreter','latex')
% title(['xlocation = ',num2str(xloc),' m'])
% niceplot_nobold_nomintick(18);
% drawnow
% pause(1)
% end 

%% extract pt along x= ct
t_itp = 0:0.1:8;
cxt_alongct= zeros(length(t_itp),length(ind_good));
ct_all = zeros(length(t_itp),length(ind_good));
for i=1:length(ind_good) %cross shore dim
    xloc = ind_good(i);
    h_xloc = h(xloc);
    c = sqrt(g*h_xloc);
    ct = c.*t_itp;
    ct_all(:,i) = ct;
    cxt_atx = squeeze(cxt(xloc,:,:));
    cxt_alongct(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct);
    
end 

cxt_alongct_mean = mean(cxt_alongct,2,'omitnan');
cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% using t_itp (s)
colormap = cmocean('phase',length(ind_good)); %parula
figure()
for i = 1:length(ind_good)
    plot(t_itp, cxt_alongct(:,i),'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
p1=plot(t_itp,cxt_alongct_mean,'LineWidth', 3,'Color','k');
hold on 
p2=plot(t_itp,cxt_alongct_25_75(:,1),'--','LineWidth', 2,'Color','k');
hold on 
p3 =plot(t_itp,cxt_alongct_25_75(:,2),'--','LineWidth', 2,'Color','k');
hold on 
clear colormap
colormap(cmocean('phase',length(ind_good)));
legend([p1, p2, p3],'Mean', '25 percentile','75 percentile')
xlabel('Time Lag (s)--along x=ct')
ylabel('Cross Correlation')
title(['Run',num2str(runnum)])
colorbar
col=colorbar;
caxis([x(ind_good(1)),x(ind_good(end))])
col.Label.String = 'Cross-shore location (m) (0 is shoreline)';
xlim([0,5])
niceplot_nobold_nomintick(18);
grid on 
hold off 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot every 10 m  
ind_subset = 1:10:length(ind_good);
colormap = colormap(cmocean('haline',length(ind_subset)));
for i = 1:length(ind_subset)
    plot(t_itp, cxt_alongct(:,ind_subset(i)),'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
clear colormap
colormap(cmocean('halin',length(ind_subset)));
xlabel('Time Lag (s)--along x=ct')
ylabel('Cross Correlation')
title(['Run',num2str(runnum)])
colorbar
col=colorbar;
caxis([x(ind_good(1)),x(ind_good(end))])
col.Label.String = 'Cross-shore location (m) (0 is shoreline)';
xlim([0,5])
niceplot_nobold_nomintick(18);
grid on 
hold off 

%% using cxt(nond x)
load('/data1/bliu/data/cxt_alongct_x.mat')
load('/data1/bliu/data/cxt_alongct_x_fitpara','fitpara','readme')

%runnum = 1;
%skipx = 5;

% ind_good = ind_good_All{runnum};
% cxt_alongct = cxt_alongct_x_ALL{runnum};
% x= x_nond_All{runnum};

runnum =1;
nond_x = fitpara(runnum).slp2.x_nond;
x_itp = 0:0.1:9;
kw = fitpara(runnum).slp2.kw;
cxt = fitpara(runnum).slp2.cxt_data;
for i  = 1:length(kw)
    nond_dx(:,i) = x_itp./kw(i);
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
xlabel('$\Delta x (m)$','Interpreter','latex')
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

%% cxt(dim x)
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

