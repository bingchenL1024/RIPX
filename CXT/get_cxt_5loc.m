% Bingchen Liu Mar 12,2025
% This code extract cxt and its fit info at 5 sz locations for paper plot 



clear
load('/data1/bliu/data/cxt_ind_good.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_withscaling.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')



ind_slp234 = [indbath.slp2;indbath.slp3;indbath.slp4];

loc5 = linspace(-0.75,-0.25,5);

for i = 1:24
    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(fitpara.slp2(i).x_nond - loc5(j))); %ind out of 68*1
    
        fitpara_5loc.slp2(i).x_nond(j) = fitpara.slp2(i).x_nond(ind_5loc_temp(j));
        fitpara_5loc.slp2(i).Tp(j) = fitpara.slp2(i).Tp(ind_5loc_temp(j));
        fitpara_5loc.slp2(i).a(j) = fitpara.slp2(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp2(i).t_itp{j} = cell2mat(fitpara.slp2(i).t_itp(ind_5loc_temp(j)));
        fitpara_5loc.slp2(i).t_nond_diva{j} = cell2mat(fitpara.slp2(i).t_itp(ind_5loc_temp(j)))./fitpara.slp2(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp2(i).cxt_data{j} = cell2mat(fitpara.slp2(i).cxt_data(ind_5loc_temp(j)));
        fitpara_5loc.slp2(i).cxt_fit{j} = cell2mat(fitpara.slp2(i).cxt_fit(ind_5loc_temp(j)));
    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(fitpara.slp3(i).x_nond - loc5(j)));
    
        fitpara_5loc.slp3(i).x_nond(j) = fitpara.slp3(i).x_nond(ind_5loc_temp(j));
        fitpara_5loc.slp3(i).Tp(j) = fitpara.slp3(i).Tp(ind_5loc_temp(j));        
        fitpara_5loc.slp3(i).a(j) = fitpara.slp3(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp3(i).t_itp{j} = cell2mat(fitpara.slp3(i).t_itp(ind_5loc_temp(j)));
        fitpara_5loc.slp3(i).t_nond_diva{j} = cell2mat(fitpara.slp3(i).t_itp(ind_5loc_temp(j)))./fitpara.slp3(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp3(i).cxt_data{j} = cell2mat(fitpara.slp3(i).cxt_data(ind_5loc_temp(j)));
        fitpara_5loc.slp3(i).cxt_fit{j} = cell2mat(fitpara.slp3(i).cxt_fit(ind_5loc_temp(j)));
    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(fitpara.slp4(i).x_nond - loc5(j)));
    
        fitpara_5loc.slp4(i).x_nond(j) = fitpara.slp4(i).x_nond(ind_5loc_temp(j));
        fitpara_5loc.slp4(i).Tp(j) = fitpara.slp4(i).Tp(ind_5loc_temp(j));        
        fitpara_5loc.slp4(i).a(j) = fitpara.slp4(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp4(i).t_itp{j} = cell2mat(fitpara.slp4(i).t_itp(ind_5loc_temp(j)));
        fitpara_5loc.slp4(i).t_nond_diva{j} = cell2mat(fitpara.slp4(i).t_itp(ind_5loc_temp(j)))./fitpara.slp4(i).a(ind_5loc_temp(j));
        fitpara_5loc.slp4(i).cxt_data{j} = cell2mat(fitpara.slp4(i).cxt_data(ind_5loc_temp(j)));
        fitpara_5loc.slp4(i).cxt_fit{j} = cell2mat(fitpara.slp4(i).cxt_fit(ind_5loc_temp(j)));
    end 



end 
%% put t at 5 loc in a vector for pdf plot 
t_decoscal_tot_5loc = []; 
for ind = 1:24
    t_decoscal_tot_5loc = [t_decoscal_tot_5loc, fitpara_5loc.slp2(ind).a,fitpara_5loc.slp3(ind).a,fitpara_5loc.slp4(ind).a];
end 
fitpara_5loc.t_deco = t_decoscal_tot_5loc;
fitpara_5loc.t_mean = mean(t_decoscal_tot_5loc,'omitmissing');
fitpara_5loc.t_std = std(t_decoscal_tot_5loc,'omitmissing');


head = 'data generated in get_cxt_5loc.m';
save('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc','fitpara_5loc','t_decoscal_tot_5loc','head')