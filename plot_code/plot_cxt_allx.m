% Bingchen Liu Aug 5, 2024
% This code plot the cxt along x= ct and it includes:
%       * all cross-shore locations
%       * every 10 m 


for runnum = 1:120
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_alongct.mat')

runpar= get_runpara(runnum);
cxt_alongct = cxt_alongct_ALL{runnum};
cxt_alongct_mean = cxt_alongct_mean{runnum};

SS = S(runnum);
x = SS.X;
x_nond = x/SS.xb;
ind_good = ind_good_All{runnum};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55 all crosshore loc
colormap = cmocean('phase',length(ind_good)); 
figure()
sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
for i = 1:length(ind_good)
    plot(t_itp, cxt_alongct(:,i),'LineWidth',1.5,'Color',colormap(i,:))
    hold on 
end 
p1=plot(t_itp,cxt_alongct_mean,'LineWidth', 3,'Color','k');
hold on 
% p2=plot(t_itp,cxt_alongct_25_75(:,1),'--','LineWidth', 2,'Color','k');
% hold on 
% p3 =plot(t_itp,cxt_alongct_25_75(:,2),'--','LineWidth', 2,'Color','k');
% hold on 
clear colormap
%colormap = cmocean('phase',length(ind_good));
colormap(cmocean('phase',length(ind_good)));

legend(p1,'Mean')
xlabel('Time Lag (s)--along x=ct')
ylabel('Cross Correlation')
colormap(cmocean('phase',length(ind_good)))
colorbar
col=colorbar;
caxis([x_nond(ind_good(1)),x_nond(ind_good(end))])
col.Label.String = 'Dimensionless cross-shore location (0 is shoreline)';
xlim([0,5])
niceplot_nobold_nomintick(18);
grid on 
hold off 


width= 15;
height = 15;
set(gcf,'Units','inches','Position',[0,0,width,height])
set(gcf,'visible','off')
saveas(gcf,['/data1/bliu/figures/RIPX_allrun/cxt/All_x/','run_',num2str(runnum),'.png'])
clf
close all
runnum
end 





