% Bingchen Liu Jan 20, 2025
% This code obtain the rmse of WB fit for all run to select the expo for WB
% fit 


%% get RMSe for 3 surfzone locations 

% clear
% rmse_diffexp=[];
% rmse_diffexp_log=[];
% rmse_diffexp_weight=[];
% rmse_diffexp_log_weight=[];
% 
% 
% badfit_num =[];
% 
% diffexpo = 1.0:0.025:1.4;
% for expo = diffexpo
% [Error_rmse] = get_rmse_WB(expo);
% rmse_diffexp = [rmse_diffexp;Error_rmse.rmse_unweight];
% rmse_diffexp_log = [rmse_diffexp_log;Error_rmse.rmse_log_unweight];
% rmse_diffexp_weight = [rmse_diffexp_weight;Error_rmse.rmse_weight];
% rmse_diffexp_log_weight = [rmse_diffexp_log_weight;Error_rmse.rmse_log_weight];
% 
% % get_weibul_qced
% % badfit_num = [badfit_num;badfit_tot];
% 
% end 
% 
% save('/data1/bliu/data/WBfit_rmseVSexpo','diffexpo','rmse_diffexp','rmse_diffexp_log','rmse_diffexp_weight','rmse_diffexp_log_weight')

%plot_rmseVSexpo

%% get RMSe for 10 surfzone locations 

clear
rmse_diffexp=[];
rmse_diffexp_log=[];
rmse_diffexp_weight=[];
rmse_diffexp_log_weight=[];


badfit_num =[];

diffexpo = 1.2:0.05:1.5;
for expo = diffexpo
tic
expo
[Error_rmse] = get_rmse_WB_10locs(expo);
rmse_diffexp = [rmse_diffexp;Error_rmse.rmse_unweight];
rmse_diffexp_log = [rmse_diffexp_log;Error_rmse.rmse_log_unweight];
rmse_diffexp_weight = [rmse_diffexp_weight;Error_rmse.rmse_weight];
rmse_diffexp_log_weight = [rmse_diffexp_log_weight;Error_rmse.rmse_log_weight];

% get_weibul_qced
% badfit_num = [badfit_num;badfit_tot];
toc
end 

save('/data1/bliu/data/WBfit_rmseVSexpo_10loc','diffexpo','rmse_diffexp','rmse_diffexp_log','rmse_diffexp_weight','rmse_diffexp_log_weight')












