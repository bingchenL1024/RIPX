% Bingchen Liu Nov 28, 2024
% This code obtain the cxt(dt) at 3 surf zone location 0.75, 0.5, 0.3 for
% plotting purpose 




clear
load('/data1/bliu/data/cxt_ind_good.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
load('/data1/bliu/data/cxt_alongct_nointerp_fitpara_qced_1mres.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')




ind_slp234 = [indbath.slp2;indbath.slp3;indbath.slp4];

loc3 = [0.75;0.5;1/3];

for i = 1:24
    for j = 1:3
        [~,ind_3loc_temp(j)] = min(abs(fitpara.slp2(i).x_nond + loc3(j)));
    
        fitpara_3loc.slp2(i).x_nond(j) = fitpara.slp2(i).x_nond(ind_3loc_temp(j));
        fitpara_3loc.slp2(i).a(j) = fitpara.slp2(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp2(i).t_itp{j} = cell2mat(fitpara.slp2(i).t_itp(ind_3loc_temp(j)));
        fitpara_3loc.slp2(i).t_nond_diva{j} = cell2mat(fitpara.slp2(i).t_itp(ind_3loc_temp(j)))./fitpara.slp2(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp2(i).cxt_data{j} = cell2mat(fitpara.slp2(i).cxt_data(ind_3loc_temp(j)));
        fitpara_3loc.slp2(i).cxt_fit{j} = cell2mat(fitpara.slp2(i).cxt_fit(ind_3loc_temp(j)));
    end 

    for j = 1:3
        [~,ind_3loc_temp(j)] = min(abs(fitpara.slp3(i).x_nond + loc3(j)));
    
        fitpara_3loc.slp3(i).x_nond(j) = fitpara.slp3(i).x_nond(ind_3loc_temp(j));
        fitpara_3loc.slp3(i).a(j) = fitpara.slp3(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp3(i).t_itp{j} = cell2mat(fitpara.slp3(i).t_itp(ind_3loc_temp(j)));
        fitpara_3loc.slp3(i).t_nond_diva{j} = cell2mat(fitpara.slp3(i).t_itp(ind_3loc_temp(j)))./fitpara.slp3(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp3(i).cxt_data{j} = cell2mat(fitpara.slp3(i).cxt_data(ind_3loc_temp(j)));
        fitpara_3loc.slp3(i).cxt_fit{j} = cell2mat(fitpara.slp3(i).cxt_fit(ind_3loc_temp(j)));
    end 

    for j = 1:3
        [~,ind_3loc_temp(j)] = min(abs(fitpara.slp4(i).x_nond + loc3(j)));
    
        fitpara_3loc.slp4(i).x_nond(j) = fitpara.slp4(i).x_nond(ind_3loc_temp(j));
        fitpara_3loc.slp4(i).a(j) = fitpara.slp4(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp4(i).t_itp{j} = cell2mat(fitpara.slp4(i).t_itp(ind_3loc_temp(j)));
        fitpara_3loc.slp4(i).t_nond_diva{j} = cell2mat(fitpara.slp4(i).t_itp(ind_3loc_temp(j)))./fitpara.slp4(i).a(ind_3loc_temp(j));
        fitpara_3loc.slp4(i).cxt_data{j} = cell2mat(fitpara.slp4(i).cxt_data(ind_3loc_temp(j)));
        fitpara_3loc.slp4(i).cxt_fit{j} = cell2mat(fitpara.slp4(i).cxt_fit(ind_3loc_temp(j)));
    end 



end 


save('/data1/bliu/data/cxt_alongct_nointerp_fitpara_qced_3loc','fitpara_3loc')