% Bingchen Liu Oct4, 2024
% This code gets the run info of parameters in CXT fit that is bad
clear

load('/data1/bliu/data/cxt_x_fitanalysis')

ind_bad = intersect(find(a_tot.slp2>1.5),find(h.slp2<1));
badrun = cell2mat(runnum_tot.slp2(ind_bad));
badrun(diff(badrun)==0) =[];


for i = 1:length(badrun)
    badruninfo(i) = get_runpara(badrun(i));
end 
