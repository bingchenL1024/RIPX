% Bingchen Liu Feb 18, 2025
% This code analyze the CXT data with dx width and for further plotting
clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/runnum_72run.mat')


%% make nond cxt vs dx width
for ind = 1:72
    run = goodrunnum(ind);
    for xloc=1:length(cxt_alongct_pm2dx_ALL{run}(1,:))
        for tlag = 1:3
            cxt_pm2dx_nond_ALL{ind}{tlag,xloc} = cxt_alongct_pm2dx_ALL{run}{tlag,xloc}./cxt_alongct_pm2dx_ALL{run}{tlag,xloc}(3); %3 corresponding to the center location
        end 
    end 
end 

%%

figure()
for runnum=1:72
    for k = 1:15:length(cxt_pm2dx_nond_ALL{goodrunnum(runnum)}(1,:))
        for i = 1:3
            plot(dxwidth_cxtcoord_pm2,cxt_pm2dx_nond_ALL{goodrunnum(runnum)}{i,k},'Color','k')
            hold on
        end 
    end 
end 
hold off
xlabel('dx lag')
ylabel('CXT/CXT(dx=0)')
ylim([-inf,inf])
niceplot_nobold(15)

figure()
for runnum=1:72
for k = 1:15:length(cxt_alongct_pm2dx_ALL{goodrunnum(runnum)}(1,:))
    for i = 2:11
        plot(dxwidth_cxtcoord_pm2,cxt_alongct_pm2dx_ALL{goodrunnum(runnum)}{i,k},'Color','k')
        hold on
    end 
end 
end 
hold off
xlabel('dx lag')
ylabel('CXT')
niceplot_nobold(15)


%%
figure()
for i = 1:72
    plot(t_itp,cxt_alongct_ALL{goodrunnum(i)},'Color','b')
    hold on 
end 
hold off