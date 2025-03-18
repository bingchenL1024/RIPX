% Bingchen Liu Mar 3, 2025
% This code analysis a b c fit from fit of cxt that didn't use
% interpolation 


clear
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_withscaling.mat')
%load('/data1/bliu/data/cxt_alongct_x_maxvar_fitpara_qced.mat')
load('/data1/bliu/data/SS_raw.mat') 


g=9.81;

% %% get index of every 5 m 
% for i = 1:120
%     ind_5mX{i,1} = intersect(S(i).X,S(i).X2);
%     ind_5mX2{i,1} = intersect(S(i).X2,S(i).X);
%     Hs_5m{i} = S(i).Hs(ind_5mX2);
% end 


%% collect data 
%%%%%%%%%%%%%%%%%%%%%%%%%% slp2
a_tot.slp2 = [];
rsq_tot.slp2 = [];
x_tot.slp2 = [];
x_br_tot.slp2 = [];
x_nond_tot.slp2 = [];
t_scale.slp2=[];
h.slp2 = [];
Tp.sl2= [];
runnum_tot.slp2 = [] ; 
Hs_interp.slp2 = [];
kw.slp2 = [];
G0_nond.slp2 = [];
for i = 1:24
    a_tot.slp2= [a_tot.slp2; fitpara.slp2(i).a];
    x_tot.slp2 = [x_tot.slp2;fitpara.slp2(i).x];
    x_br_tot.slp2 = [x_br_tot.slp2;fitpara.slp2(i).xb];
    Tp.sl2 = [Tp.sl2;fitpara.slp2(i).Tp];
    x_nond_tot.slp2 = [x_nond_tot.slp2; fitpara.slp2(i).x_nond];
    rsq_tot.slp2= [rsq_tot.slp2; fitpara.slp2(i).rsq];
    t_scale.slp2=[t_scale.slp2;fitpara.slp2(i).t_scale];
    h.slp2=[h.slp2;fitpara.slp2(i).h];
    runnum_tot.slp2 = [runnum_tot.slp2;fitpara.slp2(i).runnum] ; 
    Hs_interp.slp2 = [Hs_interp.slp2;fitpara.slp2(i).Hs_interp];
    kw.slp2 = [kw.slp2;fitpara.slp2(i).kw];
    G0_nond.slp2 = [G0_nond.slp2;fitpara.slp2(i).G0_nond];
end 

a_tot.slp2 = a_tot.slp2(find(x_nond_tot.slp2>-1));
rsq_tot.slp2 = rsq_tot.slp2(find(x_nond_tot.slp2>-1));
t_scale.slp2= t_scale.slp2(find(x_nond_tot.slp2>-1));
t_sincebr.slp2 = 2*(g*0.02)^(-0.5)*(-(-x_tot.slp2(find(x_nond_tot.slp2>-1))).^(0.5)+(x_br_tot.slp2(find(x_nond_tot.slp2>-1))).^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp3
a_tot.slp3 = [];
rsq_tot.slp3 = [];
x_tot.slp3 = [];
x_br_tot.slp3 = [];
x_nond_tot.slp3 = [];
t_scale.slp3=[];
h.slp3 = [];
Tp.sl3= [];
runnum_tot.slp3 = [] ; 
Hs_interp.slp3 = [];
kw.slp3 = [];
G0_nond.slp3 = [];

for i = 1:24
    a_tot.slp3= [a_tot.slp3; fitpara.slp3(i).a];
    % if length(a_tot.slp3) > 249
    %     keyboard
    % end 
    x_tot.slp3 = [x_tot.slp3;fitpara.slp3(i).x];
    x_br_tot.slp3 = [x_br_tot.slp3;fitpara.slp3(i).xb];
    Tp.sl3 = [Tp.sl3;fitpara.slp3(i).Tp];
    x_nond_tot.slp3 = [x_nond_tot.slp3; fitpara.slp3(i).x_nond];
    rsq_tot.slp3= [rsq_tot.slp3; fitpara.slp3(i).rsq];
    t_scale.slp3=[t_scale.slp3;fitpara.slp3(i).t_scale];
    h.slp3=[h.slp3;fitpara.slp3(i).h];
    runnum_tot.slp3 = [runnum_tot.slp3;fitpara.slp3(i).runnum] ; 
    Hs_interp.slp3 = [Hs_interp.slp3;fitpara.slp3(i).Hs_interp];
    kw.slp3 = [kw.slp3;fitpara.slp3(i).kw];
    G0_nond.slp3 = [G0_nond.slp3;fitpara.slp3(i).G0_nond];



end 

a_tot.slp3 = a_tot.slp3(find(x_nond_tot.slp3>-1));
rsq_tot.slp3 = rsq_tot.slp3(find(x_nond_tot.slp3>-1));
t_scale.slp3= t_scale.slp3(find(x_nond_tot.slp3>-1));
t_sincebr.slp3 = 2*(g*0.03)^(-0.5)*(-(-x_tot.slp3(find(x_nond_tot.slp3>-1))).^(0.5)+(x_br_tot.slp3(find(x_nond_tot.slp3>-1))).^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp4
a_tot.slp4 = [];
rsq_tot.slp4 = [];
x_tot.slp4 = [];
x_br_tot.slp4 = [];
x_nond_tot.slp4 = [];
t_scale.slp4=[];
h.slp4 = [];
Tp.sl4= [];
runnum_tot.slp4 = [];
Hs_interp.slp4 = [];
kw.slp4 = [];
G0_nond.slp4 = [];



for i = 1:24
    a_tot.slp4= [a_tot.slp4; fitpara.slp4(i).a];
    x_tot.slp4 = [x_tot.slp4;fitpara.slp4(i).x];
    x_br_tot.slp4 = [x_br_tot.slp4;fitpara.slp4(i).xb];
    Tp.sl4 = [Tp.sl4;fitpara.slp4(i).Tp];
    x_nond_tot.slp4 = [x_nond_tot.slp4; fitpara.slp4(i).x_nond];
    rsq_tot.slp4= [rsq_tot.slp4; fitpara.slp4(i).rsq];
    t_scale.slp4=[t_scale.slp4;fitpara.slp4(i).t_scale];
    h.slp4=[h.slp4;fitpara.slp4(i).h];
    runnum_tot.slp4 = [runnum_tot.slp4;fitpara.slp4(i).runnum] ; 
    Hs_interp.slp4 = [Hs_interp.slp4;fitpara.slp4(i).Hs_interp];
    kw.slp4 = [kw.slp4;fitpara.slp4(i).kw];
    G0_nond.slp4 = [G0_nond.slp4;fitpara.slp4(i).G0_nond];



end 


a_tot.slp4 = a_tot.slp4(find(x_nond_tot.slp4>-1));
rsq_tot.slp4 = rsq_tot.slp4(find(x_nond_tot.slp4>-1));
t_scale.slp4= t_scale.slp4(find(x_nond_tot.slp4>-1));
t_sincebr.slp4 = 2*(g*0.04)^(-0.5)*(-(-x_tot.slp4(find(x_nond_tot.slp4>-1))).^(0.5)+(x_br_tot.slp4(find(x_nond_tot.slp4>-1))).^0.5);




% x_nond_tot.slp2 = x_nond_tot.slp2(find(x_nond_tot.slp2>-1)); % no longer
% needed since already took care of in ind_good 
% x_nond_tot.slp3 = x_nond_tot.slp3(find(x_nond_tot.slp3>-1));
% x_nond_tot.slp4 = x_nond_tot.slp4(find(x_nond_tot.slp4>-1));

%get run index of given beach slope 
[p2,~] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,~] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,~] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

for i = 1:24
    Tp.sl2(i,1) = S(p2(i,1)).Tp_T;
    Tp.sl3(i,1) = S(p3(i,1)).Tp_T;
    Tp.sl4(i,1) = S(p4(i,1)).Tp_T;
end 


%% combine data with diff slp into one 
a_tot_all = [a_tot.slp2;a_tot.slp3;a_tot.slp4];
h_all = [h.slp2;h.slp3;h.slp4];
Tp_all = [Tp.sl2;Tp.sl3;Tp.sl4];
runnum_tot_all = cell2mat([runnum_tot.slp2;runnum_tot.slp3;runnum_tot.slp4]);
Hs_interp_tot = [Hs_interp.slp2;Hs_interp.slp3;Hs_interp.slp4];
kw_all=[kw.slp2;kw.slp3;kw.slp4];
t_scale_all = [t_scale.slp2;t_scale.slp3;t_scale.slp4];
x_nond_all = [x_nond_tot.slp2;x_nond_tot.slp3;x_nond_tot.slp4];
G0_nond_all = [G0_nond.slp2;G0_nond.slp3;G0_nond.slp4];

a_nond = a_tot_all./t_scale_all;

g=9.8;

nond_Tph = g.*Tp_all.^2./h_all;
nond_Hsh = Hs_interp_tot./h_all;
nond_Hshslp.slp2= Hs_interp.slp2./h.slp2;
nond_Hshslp.slp3= Hs_interp.slp3./h.slp3;
nond_Hshslp.slp4= Hs_interp.slp4./h.slp4;
%nond= h_all./(g.*Tp_all.^2);


%% get index for plot purpose (i.e. for diff Tp) optional
%get_runinfo_4plot_nointerp
%load('/data1/bliu/data/cxt_nointerp_runinfo_G0lim0p15.mat')
% 
% % ================================================= Tp
% for i = 1:length(runinfo_tot)
%     ind_Tp8(i,1)=  contains(runinfo_tot(i).wave,'Tp8.0');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_Tp14(i,1)=  contains(runinfo_tot(i).wave,'Tp14.0');
% end 
% 
% % ================================================= Hs
% for i = 1:length(runinfo_tot)
%     ind_Hs0p5(i,1)=  contains(runinfo_tot(i).wave,'Hs0.5');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_Hs0p8(i,1)=  contains(runinfo_tot(i).wave,'Hs0.8');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_Hs1p1(i,1)=  contains(runinfo_tot(i).wave,'Hs1.1');
% end 
% 
% % ================================================= ds
% for i = 1:length(runinfo_tot)
%     ind_ds2p5(i,1)=  contains(runinfo_tot(i).wave,'ds2.5');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_ds5(i,1)=  contains(runinfo_tot(i).wave,'ds5.0');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_ds10(i,1)=  contains(runinfo_tot(i).wave,'ds10.0');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_ds20(i,1)=  contains(runinfo_tot(i).wave,'ds20.0');
% end 
% 
% 
% % ================================================= bathy
% for i = 1:length(runinfo_tot)
%     ind_slp2(i,1)=  contains(runinfo_tot(i).bath,'002');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_slp3(i,1)=  contains(runinfo_tot(i).bath,'003');
% end 
% 
% for i = 1:length(runinfo_tot)
%     ind_slp4(i,1)=  contains(runinfo_tot(i).bath,'004');
% end 
% 
% 
% 
% 




%% Fit a linear line to get rsquare nond t vs Hs/h
ftt = strcat('A*x+B');

ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -0.5];
opts.StartPoint = [3.7 0.1]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi

[f_a_nond,gof_a_nond]=fit(Hs_interp_tot(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),...   
a_tot_all(~isnan(a_tot_all))./t_scale_all(~isnan(a_tot_all)),ft, opts);

nond_Hsh_model = 0:0.05:1;
a_nond_model = nond_Hsh_model.* f_a_nond.A+f_a_nond.B;

%% Fit nond t with nond G0


%ftt = strcat('A*exp(B*x)+C');
ftt = strcat('A*x^B+C');
%ftt = strcat('A+B*x+C*x^2');
%ftt = strcat('exp(x-A)+B');
%ftt = strcat('A*x+B');

ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -Inf];
opts.StartPoint = [1 1 0.8]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi

[f_a_nond_G0,gof_a_nond_G0]=fit(G0_nond_all(~isnan(a_tot_all)),...   
a_tot_all(~isnan(a_tot_all))./t_scale_all(~isnan(a_tot_all)),ft, opts);

nond_G0_model = 0.15:0.05:1;
%a_nond_model_G0 = f_a_nond_G0.A*exp(nond_G0_model*f_a_nond_G0.B)+f_a_nond_G0.C;
a_nond_model_G0 = f_a_nond_G0.A*nond_G0_model.^f_a_nond_G0.B+f_a_nond_G0.C;
%%
save('/data1/bliu/data/cxt_fitanalysis_max_dxwidth')

