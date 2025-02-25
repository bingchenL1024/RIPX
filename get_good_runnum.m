% Bingchen liu Feb 6, 2025
% This code get all the runs with good bathy (0.02, 0.03, 0.04)


[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

goodrunnum = [p2;p3;p4];
goodrunnum = sort(goodrunnum,'ascend');
head='run index with desired bathy(0.02,0.03,0.04). Generated from "get_good_runnum"';

save('/data1/bliu/data/runnum_72run',"head","goodrunnum")