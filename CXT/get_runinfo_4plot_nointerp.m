% Bingchen Liu Nov 13, 2024
% This code get the run parameters for cxt_nointerp data for plotting
% purpose(i.e. filter out specific run for color plot)

clear
load('/data1/bliu/data/cxt_nointerp_fitanalysis')
clear runinfo_tot


for i =1:length(runnum_tot_all)
    runinfo_tot(i)=get_runpara(runnum_tot_all(i));
end 
readme = 'created in get_runinfo_4plot_nointerp.m';
save('/data1/bliu/data/cxt_nointerp_runinfo','runinfo_tot','readme')