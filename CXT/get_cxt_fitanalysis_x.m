% Bingchen Liu Aug 15, 2024
% This code analysis a b c fit from fit of cxt 


clear
load('/data1/bliu/data/cxt_alongct_x_fitpara_qced.mat')
%load('/data1/bliu/data/cxt_alongct_x_maxvar_fitpara_qced.mat')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_runinfo')


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
b_tot.slp2 = [];
c_tot.slp2 = [];
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
for i = 1:24
    a_tot.slp2= [a_tot.slp2; fitpara.slp2(i).a];
    b_tot.slp2= [b_tot.slp2; fitpara.slp2(i).b];
    c_tot.slp2= [c_tot.slp2; fitpara.slp2(i).c];
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

end 

a_tot.slp2 = a_tot.slp2(find(x_nond_tot.slp2>-1));
b_tot.slp2 = b_tot.slp2(find(x_nond_tot.slp2>-1));
c_tot.slp2 = c_tot.slp2(find(x_nond_tot.slp2>-1));
rsq_tot.slp2 = rsq_tot.slp2(find(x_nond_tot.slp2>-1));
t_scale.slp2= t_scale.slp2(find(x_nond_tot.slp2>-1));
t_sincebr.slp2 = 2*(g*0.02)^(-0.5)*(-(-x_tot.slp2(find(x_nond_tot.slp2>-1))).^(0.5)+(x_br_tot.slp2(find(x_nond_tot.slp2>-1))).^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp3
a_tot.slp3 = [];
b_tot.slp3 = [];
c_tot.slp3 = [];
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

for i = 1:24
    a_tot.slp3= [a_tot.slp3; fitpara.slp3(i).a];
    b_tot.slp3= [b_tot.slp3; fitpara.slp3(i).b];
    c_tot.slp3= [c_tot.slp3; fitpara.slp3(i).c];
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


end 

a_tot.slp3 = a_tot.slp3(find(x_nond_tot.slp3>-1));
b_tot.slp3 = b_tot.slp3(find(x_nond_tot.slp3>-1));
c_tot.slp3 = c_tot.slp3(find(x_nond_tot.slp3>-1));
rsq_tot.slp3 = rsq_tot.slp3(find(x_nond_tot.slp3>-1));
t_scale.slp3= t_scale.slp3(find(x_nond_tot.slp3>-1));
t_sincebr.slp3 = 2*(g*0.03)^(-0.5)*(-(-x_tot.slp3(find(x_nond_tot.slp3>-1))).^(0.5)+(x_br_tot.slp3(find(x_nond_tot.slp3>-1))).^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp4
a_tot.slp4 = [];
b_tot.slp4 = [];
c_tot.slp4 = [];
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



for i = 1:24
    a_tot.slp4= [a_tot.slp4; fitpara.slp4(i).a];
    b_tot.slp4= [b_tot.slp4; fitpara.slp4(i).b];
    c_tot.slp4= [c_tot.slp4; fitpara.slp4(i).c];
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



end 


a_tot.slp4 = a_tot.slp4(find(x_nond_tot.slp4>-1));
b_tot.slp4 = b_tot.slp4(find(x_nond_tot.slp4>-1));
c_tot.slp4 = c_tot.slp4(find(x_nond_tot.slp4>-1));
rsq_tot.slp4 = rsq_tot.slp4(find(x_nond_tot.slp4>-1));
t_scale.slp4= t_scale.slp4(find(x_nond_tot.slp4>-1));
t_sincebr.slp4 = 2*(g*0.04)^(-0.5)*(-(-x_tot.slp4(find(x_nond_tot.slp4>-1))).^(0.5)+(x_br_tot.slp4(find(x_nond_tot.slp4>-1))).^0.5);




x_nond_tot.slp2 = x_nond_tot.slp2(find(x_nond_tot.slp2>-1));
x_nond_tot.slp3 = x_nond_tot.slp3(find(x_nond_tot.slp3>-1));
x_nond_tot.slp4 = x_nond_tot.slp4(find(x_nond_tot.slp4>-1));

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
b_tot_all = [b_tot.slp2;b_tot.slp3;b_tot.slp4];
c_tot_all = [c_tot.slp2;c_tot.slp3;c_tot.slp4];
h_all = [h.slp2;h.slp3;h.slp4];
Tp_all = [Tp.sl2;Tp.sl3;Tp.sl4];
runnum_tot_all = cell2mat([runnum_tot.slp2;runnum_tot.slp3;runnum_tot.slp4]);
Hs_interp_tot = [Hs_interp.slp2;Hs_interp.slp3;Hs_interp.slp4];
kw_all=[kw.slp2;kw.slp3;kw.slp4];

g=9.8;

nond = g.*Tp_all.^2./h_all;
nond2 = Hs_interp_tot./h_all;
%nond= h_all./(g.*Tp_all.^2);
% ================================================= Tp
for i = 1:length(runinfo_tot)
    ind_Tp8(i,1)=  contains(runinfo_tot(i).wave,'Tp8.0');
end 

for i = 1:length(runinfo_tot)
    ind_Tp14(i,1)=  contains(runinfo_tot(i).wave,'Tp14.0');
end 

% ================================================= Hs
for i = 1:length(runinfo_tot)
    ind_Hs0p5(i,1)=  contains(runinfo_tot(i).wave,'Hs0.5');
end 

for i = 1:length(runinfo_tot)
    ind_Hs0p8(i,1)=  contains(runinfo_tot(i).wave,'Hs0.8');
end 

for i = 1:length(runinfo_tot)
    ind_Hs1p1(i,1)=  contains(runinfo_tot(i).wave,'Hs1.1');
end 

% ================================================= ds
for i = 1:length(runinfo_tot)
    ind_ds2p5(i,1)=  contains(runinfo_tot(i).wave,'ds2.5');
end 

for i = 1:length(runinfo_tot)
    ind_ds5(i,1)=  contains(runinfo_tot(i).wave,'ds5.0');
end 

for i = 1:length(runinfo_tot)
    ind_ds10(i,1)=  contains(runinfo_tot(i).wave,'ds10.0');
end 

for i = 1:length(runinfo_tot)
    ind_ds20(i,1)=  contains(runinfo_tot(i).wave,'ds20.0');
end 


% ================================================= bathy
for i = 1:length(runinfo_tot)
    ind_slp2(i,1)=  contains(runinfo_tot(i).bath,'002');
end 

for i = 1:length(runinfo_tot)
    ind_slp3(i,1)=  contains(runinfo_tot(i).bath,'003');
end 

for i = 1:length(runinfo_tot)
    ind_slp4(i,1)=  contains(runinfo_tot(i).bath,'004');
end 


%% Fit a(b) VS h
ftt = strcat('A*x+B');

ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -0.5];
opts.StartPoint = [0.5 0.5]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi
[f_a,gof_a]=fit(h_all(~isnan(a_tot_all)),a_tot_all(~isnan(a_tot_all)),ft, opts);
[f_b,gof_b]=fit(h_all(~isnan(b_tot_all)),b_tot_all(~isnan(b_tot_all)),ft, opts);

opts.StartPoint = [1.8 0.5]; % beginning parameters - amp, mu, std.
[f_a_nond,gof_a_nond]=fit(Hs_interp_tot(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),...   
a_tot_all(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),ft, opts);

opts.StartPoint = [1 0.25]; % beginning parameters - amp, mu, std.
[f_b_nond,gof_b_nond]=fit(Hs_interp_tot(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),...   
b_tot_all(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),ft, opts);

opts.StartPoint = [-1 1]; 
[f_c_nond,gof_c_nond]=fit(Hs_interp_tot(~isnan(a_tot_all))./h_all(~isnan(a_tot_all)),...   
c_tot_all(~isnan(a_tot_all)),ft, opts);

h_model = 0:0.1:4;
a_model = h_model.*f_a.A+f_a.B;
b_model = h_model.*f_b.A+f_b.B;

Hs_nond_model = 0:0.1:1.5;
a_nond_model = Hs_nond_model.* f_a_nond.A+f_a_nond.B;
b_nond_model = Hs_nond_model.* f_b_nond.A+f_b_nond.B;
c_nond_model = Hs_nond_model.* f_c_nond.A+f_c_nond.B;


%%
save('/data1/bliu/data/cxt_x_fitanalysis')
%save('/data1/bliu/data/cxt_x_fitanalysis','a_tot','b_tot','c_tot','x_nond_tot','x_tot','t_scale','h','Tp','runnum_tot','Hs_interp')
%save('/data1/bliu/data/cxt_x_maxvar_fitanalysis','a_tot','b_tot','c_tot','x_nond_tot','x_tot','t_scale','h','Tp','runnum_tot')

%% relation between a and b 

% include = find(a_tot.slp2>0.1 & a_tot.slp2<1&b_tot.slp2<0.5&b_tot.slp2>0);
% 
% %ftt = strcat('real(acos(a*exp(1/(x))))') ;
% ftt = strcat('real(b*acos(a*exp(1/(x)))+c)') ;
% ft = fittype( sprintf('%s',ftt));
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.Lower = [-20 -20 -20];
% opts.StartPoint = [-0.1 0 0]; % beginning parameters - amp, mu, std.
% opts.Upper = [0 20 20];
% [f,gof]=fit(a_tot.slp2(include),b_tot.slp2(include),ft, opts);
% a_fit.slp2= 0:0.1:1.5;
% b_fit.slp2 = real(f.b.*acos(f.a.*exp(1./(a_fit.slp2)))+f.c);
% %b_fit = real(acos(f.a.*exp(1./a_fit)));
% 
% 
% ftt1 = strcat('1/(c*(x+a))+b');
% ft1 = fittype( sprintf('%s',ftt1));
% opts1 = fitoptions( ft1 );
% opts1.Display = 'Off';
% opts1.Lower = [-inf -20 -inf];
% opts1.StartPoint = [-0.1 0 1]; % beginning parameters - amp, mu, std.
% opts1.Upper = [-20 20 inf];
% [f1,gof1]=fit(a_tot.slp2(include),b_tot.slp2(include),ft1, opts1);
% a_fit.slp2= 0:0.1:1.5;
% b_fit1.slp2 = (1./(f1.c.*(a_fit.slp2+f1.a))) + f1.b;
% 
% figure()
% scatter(a_tot.slp2(include),b_tot.slp2(include),50,'filled','b')
% hold on 
% plot(a_fit.slp2,b_fit.slp2,'LineWidth',2)
% %plot(a_fit.slp2,b_fit1.slp2,'LineWidth',2)
% xlabel('a')
% ylabel('b')
% niceplot_nobold_nomintick(18)
% legend('data','fit')
% title('Using fit = real(b*acos(a*exp(1/(x)))+c)')
% hold off 
% %ylim([0 0.5])
% 

