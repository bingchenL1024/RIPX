% Bingchen Liu Mar 14, 2025
% This code compute the scaling for the fit parameter \tau and add it to
% variable 'fitpara'
% July 4: add cFbr_nond scaling test 

clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_1mres.mat')

load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/ind_of_diff_bath.mat')
load('/data1/bliu/data/cxt_ind_good.mat')


xskip = 1; %%NOTE!!!!: it needs to be consistent with xskip in get_cxt_fit_nointerp_max_dxwidth
g=9.81;
t_scaleexpr = '(h/g).^(1/2)';
%t_scaleexpr2 = 'g.*h.^2.*Dw_interp.^(-1)';
%t_scaleexpr2 = 'h.*Dw_interp.^(-1/3)';
t_scaleexpr2 = 'Dw_interp.^(1/3)/g';
%t_scaleexpr2 = 'dirspr_loc';
cFbr_symbol = 'Dw_interp./((g*h).^(0.5)*SS.hb^2)';


for ind_slp2 = 1:24
    runnum = indbath.slp2(ind_slp2);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    
    beta= 0.02;
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    Ir_br= 0.02/(SS.stpb)^0.5;
    Dw_interp = real(interp1(SS.X2,SS.dECG,x)); 
    dirspr_loc = interp1(SS.X2,SS.sigma_th,x)/SS.sigma_b;


    %t_scale = sqrt(h/g);
    t_scale = eval(t_scaleexpr);
    t_scale2 = eval(t_scaleexpr2);

    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification
    curlFbr_para = eval(cFbr_symbol);
    curlFbr_para_nond = curlFbr_para/max(curlFbr_para);

    j = 0;
    for i = 1:xskip:dim_cxt(2) %crossshore scale ----------fit every 1 meter
        j = j+1;      
        fitpara.slp2(ind_slp2,1).t_scale(j,1) = t_scale(i);  
        fitpara.slp2(ind_slp2,1).t_scale2(j,1) = t_scale2(i);  
        fitpara.slp2(ind_slp2,1).curlFbr_para_nond(j,1) = curlFbr_para_nond(i); 
        fitpara.slp2(ind_slp2,1).dirspr_loc(j,1) = dirspr_loc(i); 

    end 
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp3

for ind_slp3 = 1:24
    runnum = indbath.slp3(ind_slp3);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);

    beta= 0.03;
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    Ir_br= 0.03/(SS.stpb)^0.5;
    Dw_interp = real(interp1(SS.X2,SS.dECG,x)); 
    dirspr_loc = interp1(SS.X2,SS.sigma_th,x)/SS.sigma_b;

    %t_scale = sqrt(h/g);
    t_scale = eval(t_scaleexpr);
    t_scale2 = eval(t_scaleexpr2);


    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification
    curlFbr_para = eval(cFbr_symbol);
    curlFbr_para_nond = curlFbr_para/max(curlFbr_para);

    j = 0;
    for i = 1:xskip:dim_cxt(2) 
        j = j+1;      
        fitpara.slp3(ind_slp3,1).t_scale(j,1) = t_scale(i); 
        fitpara.slp3(ind_slp3,1).t_scale2(j,1) = t_scale2(i);  
        fitpara.slp3(ind_slp3,1).curlFbr_para_nond(j,1) = curlFbr_para_nond(i); 
        fitpara.slp3(ind_slp3,1).dirspr_loc(j,1) = dirspr_loc(i); 

    end 
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp4
for ind_slp4 = 1:24
    runnum = indbath.slp4(ind_slp4);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    Ir_br= 0.04/(SS.stpb)^0.5;
    Dw_interp = real(interp1(SS.X2,SS.dECG,x)); 
    dirspr_loc = interp1(SS.X2,SS.sigma_th,x)/SS.sigma_b;

    beta = 0.04;
    %t_scale = sqrt(h/g);
    t_scale = eval(t_scaleexpr);
    t_scale2 = eval(t_scaleexpr2);



    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification
    curlFbr_para = eval(cFbr_symbol);
    curlFbr_para_nond = curlFbr_para/max(curlFbr_para);

    j = 0;
    for i = 1:xskip:dim_cxt(2) 
        j = j+1;            
        fitpara.slp4(ind_slp4,1).t_scale(j,1) = t_scale(i);
        fitpara.slp4(ind_slp4,1).t_scale2(j,1) = t_scale2(i);  
        fitpara.slp4(ind_slp4,1).curlFbr_para_nond(j,1) = curlFbr_para_nond(i); 
        fitpara.slp4(ind_slp4,1).dirspr_loc(j,1) = dirspr_loc(i); 

    end 
end



head = 'generated in get_CXT_scaling2fit';
save('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_withscaling.mat','fitpara')
