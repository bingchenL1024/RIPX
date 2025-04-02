% Bingchen Liu Mar 11, 2025 
% This code extract 4 CXT at specific x for 4 panels plot 


clear
close all

load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen_qced.mat')
load('/data1/bliu/data/runnum_72run')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
load('/data1/bliu/data/ind_of_diff_bath.mat')
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_withscaling.mat')

%%
runnum=17;
x_nond_ind = linspace(-0.625,-0.25,4);


g = 9.81;



dim = size(cell2mat(CXT_ALL(runnum)));
SS = S(runnum);
x = SS.X;
x2= SS.X2;
h = SS.h;
Hs = SS.Hs;
Hs_interp = interp1(x2,Hs,x);
cxt = cell2mat(CXT_ALL(runnum));
cxt_upbd=cell2mat(CXT_upbd_ALL(runnum));
cxt_lowbd=cell2mat(CXT_lowbd_ALL(runnum));
xb= SS.xb;
x_nond= x_nond_All{runnum};
%x_nond = x./xb;
dx_max = dx_max_All{runnum};
G0_nond = abs(SS.curlF_std)./abs(max(SS.curlF_std));%G0 modification
ind_good = ind_good_All{runnum};
c_fit_1run = c_fit{runnum};

for j =1:4 % 4 cross-shore locations
[~,ind_4loc(j)] = min(abs(x_nond - x_nond_ind(j)));
%[~,ind_4loc_fit(j)] = min(abs(fitpara.slp2(i).x_nond - x_nond_ind(j))); %ind out of 68*1

i=ind_4loc(j);%index for ind_good (68*1) vector
xind = ind_good(i); %ind for allx (147*1) vector convert it to all x coordinate
dx_max_loc(:,j) = dx_max(:,i);
dt_max_loc = 0:5;
h_xloc = h(xind);
Hs_xloc= Hs_interp(xind);
c(j) = sqrt(g*h_xloc);
c_fit_4loc(j) = c_fit_1run(i);
cxt_atx(:,:,j) = squeeze(cxt(xind,:,:));

end 

head = 'data generated in get_4panel_cxtatx.m';

save('/data1/bliu/data/cxt_atx_4panels','cxt_atx','c','c_fit_4loc','dt_max_loc','dx_max_loc','x_nond_ind','head')
