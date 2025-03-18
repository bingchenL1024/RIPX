% Bingchen Liu Mar 3, 2025
% This code calculate the nond dx width and to for plots of CXT vs dx width
% and nond CXT vs nond dx width

%dx_nond1 = dx/(c \tao(deco t scale))
%dx_nond2 = dx/(c dt(time lag))

clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/runnum_72run.mat')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')

dx_width = -2:1:2;

for ind = 1:24
    runnum = indbath.slp2(ind);
    a_decotscale = fitpara.slp2(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp2{ind,1}.dx_nond1{1,xloc} = dx_width/(c_tot{runnum}(xloc)*fitpara.slp2(ind).a(xloc));
        cxt_dxwidth.slp2{ind,1}.dx_nond2{1,xloc} = dx_width/(c_tot{runnum}(xloc));
        cxt_dxwidth.slp2{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp2{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

    end 
end 


for ind = 1:24
    runnum = indbath.slp3(ind);
    a_decotscale = fitpara.slp3(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp3{ind,1}.dx_nond1{1,xloc} = dx_width/(c_tot{runnum}(xloc)*fitpara.slp3(ind).a(xloc));
        cxt_dxwidth.slp3{ind,1}.dx_nond2{1,xloc} = dx_width/(c_tot{runnum}(xloc));
        cxt_dxwidth.slp3{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp3{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

    end 
end 

for ind = 1:24
    runnum = indbath.slp4(ind);
    a_decotscale = fitpara.slp4(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp4{ind,1}.dx_nond1{1,xloc} = dx_width/(c_tot{runnum}(xloc)*fitpara.slp4(ind).a(xloc));
        cxt_dxwidth.slp4{ind,1}.dx_nond2{1,xloc} = dx_width/(c_tot{runnum}(xloc)*fitpara.slp4(ind).a(xloc));
        cxt_dxwidth.slp4{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp4{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

    end 
end 


cxt_dxwidth.slp2  = cell2mat(cxt_dxwidth.slp2);
cxt_dxwidth.slp3  = cell2mat(cxt_dxwidth.slp3);
cxt_dxwidth.slp4  = cell2mat(cxt_dxwidth.slp4);


save('/data1/bliu/data/cxt_dxwidth',"cxt_dxwidth")


