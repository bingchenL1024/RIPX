% Bingchen Liu June 20, 2025
% this code explore different scaling and its effect on final data collapse
% 
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


%% final scaling 
% ========================== take the good one for final plot

g = 9.81;

kh.slp2 = SS.h2.*SS.kl2pm;
kh.slp3 = SS.h3.*SS.kl3pm;
kh.slp4 = SS.h4.*SS.kl4pm;

Sk.Sk2_hb_Ir_div = (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5))).*SS.Irb2;
Sk.Sk3_hb_Ir_div = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5))).*SS.Irb3;
Sk.Sk4_hb_Ir_div = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5))).*SS.Irb4;

Sk.Sk2_hb_beta= (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5))).*0.02;
Sk.Sk3_hb_beta = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5))).*0.03;
Sk.Sk4_hb_beta = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5))).*0.04;

Sk.Sk2_hb = (SS.Dw2./(SS.hb2.^2.*((g.*(SS.h2)).^0.5)));
Sk.Sk3_hb = (SS.Dw3./(SS.hb3.^2.*((g.*(SS.h3)).^0.5)));
Sk.Sk4_hb = (SS.Dw4./(SS.hb4.^2.*((g.*(SS.h4)).^0.5)));


G0_nond_Ir.slp2 = sqrt(A.wfit2.intw)./Sk.Sk2_hb_Ir_div;
G0_nond_Ir.slp3 = sqrt(A.wfit3.intw)./Sk.Sk3_hb_Ir_div;
G0_nond_Ir.slp4 = sqrt(A.wfit4.intw)./Sk.Sk4_hb_Ir_div;

G0_nond_beta.slp2 = sqrt(A.wfit2.intw)./Sk.Sk2_hb_beta;
G0_nond_beta.slp3 = sqrt(A.wfit3.intw)./Sk.Sk3_hb_beta;
G0_nond_beta.slp4 = sqrt(A.wfit4.intw)./Sk.Sk4_hb_beta;

G0_nond.slp2 = sqrt(A.wfit2.intw)./Sk.Sk2_hb;
G0_nond.slp3 = sqrt(A.wfit3.intw)./Sk.Sk3_hb;
G0_nond.slp4 = sqrt(A.wfit4.intw)./Sk.Sk4_hb;

%% calculate rsquare 
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');

loc_num=length(G0.G0_nond.slp2(1,:));
x_num = length(G0_nond.slp2(1,:));

ds_10loc.slp2 =repmat(G0.ds.sigtb2,[loc_num,1])';
ds_10loc.slp3 =repmat(G0.ds.sigtb3,[loc_num,1])';
ds_10loc.slp4 =repmat(G0.ds.sigtb4,[loc_num,1])';

rsq = get_rsq(G0_nond,ds_10loc)
rsq_beta = get_rsq(G0_nond_beta,ds_10loc)
rsq_Ir = get_rsq(G0_nond_Ir,ds_10loc)


%% plot 
markersz=35;


col = cmocean('thermal',x_num);


figure()
subplot(131)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0_nond.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),G0_nond.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),G0_nond.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$\frac{\hat{G}_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} $','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)

subplot(132)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0_nond_beta.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),G0_nond_beta.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),G0_nond_beta.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$\frac{\hat{G}_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} \beta^{-1}$','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)

subplot(133)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),G0_nond_Ir.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),G0_nond_Ir.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),G0_nond_Ir.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$\frac{\hat{G}_{0} (gh)^{1/2} h_{b}^2}{D_\mathrm{w}} Ir^{-1}$','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)




%% function 

function rsq = get_rsq(G0_nond,ds)

G0_nond_all = [G0_nond.slp2(:);G0_nond.slp3(:);G0_nond.slp4(:);];
ds_all = [ds.slp2(:),ds.slp3(:),ds.slp4(:)];
ftt = strcat('A*x+B');
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -0.5];
opts.StartPoint = [0.5 2]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi

[~,gof]=fit(ds_all(~isnan(G0_nond_all)),G0_nond_all(~isnan(G0_nond_all)),ft, opts);
rsq=gof.rsquare;
end 