% Bingchen Liu May 20, 2024
% This code calculated all the variables/data needed for the
% "plot_G0_scaling_3loc.m" code
% This code calculate the scale for G0 and save the scales for nondF plot
% This code also include data 
% "K": include directional spread and peak wave number for incoming wave
% spectra
% "A": spectra info from model data and Weibull fit
% "SS": basic model run info from S.mat
% This analysis code only process 3 surf zone locations 


% note '_2025' indicates this is post AGU2024 version and save for paper
% plot in different file while keep the AGU2024 figure data

clear

%addpath(genpath('/home/ssuanda/matlab/')) %modified by BL 
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
addpath(genpath('/data1/bliu'))
addpath(genpath('/home/ffeddersen/matlab'))


%load('/data1/bliu/data/qa_qc_RIPX_NK.mat')
load('/data1/bliu/data/SS_A_ds_qced_10loc')

%load qa_qc_RIPX_NK.mat
%load edited_Wall.mat   %original code 
%SS = load('/data1/bliu/data/S0scale_10');

loc_num = length(A.wfit2.intw(1,:));

%% 
%for plotting
% fs = 12;
% col2 = [0.3 0.3 1;0.9 0.6 0;1 .2 .2];
% 
% for plotting surf zone cross shore location 
% sz_loc = [-0.75, -0.5, -0.33];
% 
% 
% beta = 1.2; % Weibull fit exponent

%% try different scaling 
% phase velocity of solitary wave
% c_solit2 = sqrt(g.*(SS.h2+SS.Hs2));
% c_solit3 = sqrt(g.*(SS.h3+SS.Hs3));
% c_solit4 = sqrt(g.*(SS.h4+SS.Hs4));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nond number
% 
% gamma.slp2= SS.Hs2./SS.h2;
% gamma.slp3= SS.Hs3./SS.h3;
% gamma.slp4= SS.Hs4./SS.h4;
% 
% 
% % best so far 
% Sk.Sk2_hb_Ir_div_solit = (SS.Dw2./(SS.hb2.^2.*(c_solit2)))./SS.Ir2;
% Sk.Sk3_hb_Ir_div_solit = (SS.Dw3./(SS.hb3.^2.*(c_solit3)))./SS.Ir3;
% Sk.Sk4_hb_Ir_div_solit = (SS.Dw4./(SS.hb4.^2.*(c_solit4)))./SS.Ir4;
% 
% Sk.Sk2_hb_Ir_div_solit_gamma = (SS.Dw2./(SS.hb2.^2.*(c_solit2)))./(SS.Ir2.*gamma.slp2.^3);
% Sk.Sk3_hb_Ir_div_solit_gamma = (SS.Dw3./(SS.hb3.^2.*(c_solit3)))./(SS.Ir3.*gamma.slp3.^3);
% Sk.Sk4_hb_Ir_div_solit_gamma = (SS.Dw4./(SS.hb4.^2.*(c_solit4)))./(SS.Ir4.*gamma.slp4.^3);
% 
% % Try different scales:
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% using solitary wave phase velocity
% 
% Sk.Sk2_c_solit = SS.Dw2./(SS.hb2.*SS.h2.*c_solit2);
% Sk.Sk3_c_solit = SS.Dw3./(SS.hb3.*SS.h3.*c_solit3);
% Sk.Sk4_c_solit = SS.Dw4./(SS.hb4.*SS.h4.*c_solit4);
% 
% Sk.Sk2_c_solit_beta_div = (SS.Dw2./(SS.hb2.*SS.h2.*c_solit2))./SS.beta2;
% Sk.Sk3_c_solit_beta_div = (SS.Dw3./(SS.hb3.*SS.h3.*c_solit3))./SS.beta3;
% Sk.Sk4_c_solit_beta_div = (SS.Dw4./(SS.hb4.*SS.h4.*c_solit4))./SS.beta4;
% 
% Sk.Sk2_c_solit_Ir_div = (SS.Dw2./(SS.hb2.*SS.h2.*c_solit2))./SS.Ir2;
% Sk.Sk3_c_solit_Ir_div = (SS.Dw3./(SS.hb3.*SS.h3.*c_solit3))./SS.Ir3;
% Sk.Sk4_c_solit_Ir_div = (SS.Dw4./(SS.hb4.*SS.h4.*c_solit4))./SS.Ir4;
% 
% Sk.Sk2_hb_Ir_div_solit = (SS.Dw2./(SS.hb2.^2.*(c_solit2)))./SS.Ir2;
% Sk.Sk3_hb_Ir_div_solit = (SS.Dw3./(SS.hb3.^2.*(c_solit3)))./SS.Ir3;
% Sk.Sk4_hb_Ir_div_solit = (SS.Dw4./(SS.hb4.^2.*(c_solit4)))./SS.Ir4;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% using hb
% Sk.Sk2_hb = SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5));
% Sk.Sk3_hb = SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5));
% Sk.Sk4_hb = SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5));
% 
% Sk.Sk2_hb_beta_div = (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5)))./SS.beta2;
% Sk.Sk3_hb_beta_div = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5)))./SS.beta3;
% Sk.Sk4_hb_beta_div = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5)))./SS.beta4;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using h
% Sk.Sk2_h = SS.Dw2./(SS.h2.*((g.*(SS.h2.^3)).^.5));
% Sk.Sk3_h = SS.Dw3./(SS.h3.*((g.*(SS.h3.^3)).^.5));
% Sk.Sk4_h = SS.Dw4./(SS.h4.*((g.*(SS.h4.^3)).^.5));
% 
% Sk.Sk2_h_Ir_div = (SS.Dw2./((SS.h2.*((g.*(SS.h2.^3)).^.5))))./SS.Ir2; %include Iribarren number
% Sk.Sk3_h_Ir_div = (SS.Dw3./((SS.h3.*((g.*(SS.h3.^3)).^.5))))./SS.Ir3; %include Iribarren number
% Sk.Sk4_h_Ir_div = (SS.Dw4./((SS.h4.*((g.*(SS.h4.^3)).^.5))))./SS.Ir4; %include Iribarren number
% 
% Sk.Sk2_h_Ir_X = (SS.Dw2./((SS.h2.*((g.*(SS.h2.^3)).^.5)))).*SS.Ir2; %include Iribarren number
% Sk.Sk3_h_Ir_X = (SS.Dw3./((SS.h3.*((g.*(SS.h3.^3)).^.5)))).*SS.Ir3; %include Iribarren number
% Sk.Sk4_h_Ir_X = (SS.Dw4./((SS.h4.*((g.*(SS.h4.^3)).^.5)))).*SS.Ir4; %include Iribarren number
% 
% 
% Sk.Sk2_h_beta_div = (SS.Dw2./((SS.h2.*((g.*(SS.h2.^3)).^.5))))./SS.beta2; %include Iribarren number
% Sk.Sk3_h_beta_div = (SS.Dw3./((SS.h3.*((g.*(SS.h3.^3)).^.5))))./SS.beta3; %include Iribarren number
% Sk.Sk4_h_beta_div = (SS.Dw4./((SS.h4.*((g.*(SS.h4.^3)).^.5))))./SS.beta4; %include Iribarren number
% 
% 
% Sk.Sk2_h_beta_X = (SS.Dw2./((SS.h2.*((g.*(SS.h2.^3)).^.5)))).*SS.beta2; %include Iribarren number
% Sk.Sk3_h_beta_X = (SS.Dw3./((SS.h3.*((g.*(SS.h3.^3)).^.5)))).*SS.beta3; %include Iribarren number
% Sk.Sk4_h_beta_X = (SS.Dw4./((SS.h4.*((g.*(SS.h4.^3)).^.5)))).*SS.beta4; %include Iribarren number
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using kw
% Sk.Sk2_k = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5));
% Sk.Sk3_k = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5));
% Sk.Sk4_k = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5));
% 
% Sk.Sk2_k_Ir_div = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5))./SS.Ir2;
% Sk.Sk3_k_Ir_div = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5))./SS.Ir3;
% Sk.Sk4_k_Ir_div = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5))./SS.Ir4;
% 
% Sk.Sk2_k_Ir_X = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5)).*SS.Ir2;
% Sk.Sk3_k_Ir_X = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5)).*SS.Ir3;
% Sk.Sk4_k_Ir_X = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5)).*SS.Ir4;
% 
% Sk.Sk2_k_beta_div = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5))./SS.beta2;
% Sk.Sk3_k_beta_div = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5))./SS.beta3;
% Sk.Sk4_k_beta_div = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5))./SS.beta4;
% 
% Sk.Sk2_k_beta_X = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5)).*SS.beta2;
% Sk.Sk3_k_beta_X = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5)).*SS.beta3;
% Sk.Sk4_k_beta_X = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5)).*SS.beta4;
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adding kh to the scale 
% 
% 
% % ============================ with kw 
% Sk.Sk2_k_kh = (SS.kl2pm.*SS.Dw2./((g.*(SS.h2.^3)).^.5)).*(kh.slp2).^4;
% Sk.Sk3_k_kh = (SS.kl3pm.*SS.Dw3./((g.*(SS.h3.^3)).^.5)).*(kh.slp3).^4;
% Sk.Sk4_k_kh = (SS.kl4pm.*SS.Dw4./((g.*(SS.h4.^3)).^.5)).*(kh.slp4).^4;
% 
% % ============================ with h
% Sk.Sk2_h_kh = SS.Dw2./(SS.h2.*((g.*(SS.h2.^3)).^5));
% Sk.Sk3_h_kh = SS.Dw3./(SS.h3.*((g.*(SS.h3.^3)).^5));
% Sk.Sk4_h_kh = SS.Dw4./(SS.h4.*((g.*(SS.h4.^3)).^5));
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adding gamma to the scale 
% Sk.Sk2_gamma_k = (SS.Dw2.*SS.kl2pm)./(((g.*SS.h2.^3).^0.5).*gamma.slp2.^3);
% Sk.Sk3_gamma_k = (SS.Dw3.*SS.kl3pm)./(((g.*SS.h3.^3).^0.5).*gamma.slp3.^3);
% Sk.Sk4_gamma_k = (SS.Dw4.*SS.kl4pm)./(((g.*SS.h4.^3).^0.5).*gamma.slp4.^3);

%% final scaling 
% ========================== take the good one for final plot

g = 9.81;

kh.slp2 = SS.h2.*SS.kl2pm;
kh.slp3 = SS.h3.*SS.kl3pm;
kh.slp4 = SS.h4.*SS.kl4pm;

% Sk.Sk2_hb_Ir_div = (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5)))./SS.Irb2;
% Sk.Sk3_hb_Ir_div = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5)))./SS.Irb3;
% Sk.Sk4_hb_Ir_div = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5)))./SS.Irb4;

Sk.Sk2_hb_beta = (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5)))./0.02;
Sk.Sk3_hb_beta = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5)))./0.03;
Sk.Sk4_hb_beta = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5)))./0.04;

G0_nond.slp2 = sqrt(A.wfit2.intw)./Sk.Sk2_hb_beta;
G0_nond.slp3 = sqrt(A.wfit3.intw)./Sk.Sk3_hb_beta;
G0_nond.slp4 = sqrt(A.wfit4.intw)./Sk.Sk4_hb_beta;

G0_data.slp2 = sqrt(A.wfit2.intw);
G0_data.slp3 = sqrt(A.wfit3.intw);
G0_data.slp4 = sqrt(A.wfit4.intw);% G0_nond.slp2 = sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div_solit;
% G0_nond.slp3 = sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div_solit;
% G0_nond.slp4 = sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div_solit;
% 



G0_nond_all = [G0_nond.slp2(:);G0_nond.slp3(:);G0_nond.slp4(:);];
ds_all = [repmat(ds.sigtb2',loc_num,1);repmat(ds.sigtb3',loc_num,1);repmat(ds.sigtb4',loc_num,1)];


ftt = strcat('A*x+B');
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -0.5];
opts.StartPoint = [0.5 2]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi

[f,gof]=fit(ds_all(~isnan(G0_nond_all)),G0_nond_all(~isnan(G0_nond_all)),ft, opts);

ds_model = 0:0.1:15;
G0_model = ds_model.* f.A+f.B;
%% save 

README = "Data generated from 'get_plotready_G0nond_10loc.m'";
README = [README "Plot ready data, including AFTER qc 'ed'SS','A','ds','Sk'(scaling)"];
README = splitlines(README)';



save('/data1/bliu/data/plotready_G0_nond_10loc_2025','README','G0_nond','G0_nond_all','G0_model','ds','ds_all','ds_model','f','gof','G0_data')
save('/data1/bliu/data/G0_nond_10loc_allsk_2025','README','SS','ds','A','Sk','kh')