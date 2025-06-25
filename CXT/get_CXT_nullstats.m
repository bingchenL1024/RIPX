% Bingchen Liu June 20, 2025
% This code obtains the stats of CXT at supposed 0 locations to determine
% the cut off


clear
load('/data1/bliu/data/runnum_72run.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')

dt_zero= 6:8; %tlag 3s-6s
dx_zero= 2:4; % xlag 0m-3m

cxt_null_tot = [];


for ind = 1:72
    runnum = goodrunnum(ind);
    cxt_1run = CXT_ALL{runnum};
    dim = size(cxt_1run);
    for xind = 1:dim(1)
        cxt_snap= squeeze(cxt_1run(xind,:,:));
        cxt_null_snap = cxt_snap(dx_zero,dt_zero);
        cxt_null_tot = [cxt_null_tot;cxt_null_snap(:)];
    end 
end 

cxt_null_tot = cxt_null_tot(~isnan(cxt_null_tot));
cxt_null_tot = sort(abs(cxt_null_tot),'descend');


threshold = 0.99;
ind_threshold = ceil(length(cxt_null_tot)*(1-threshold));
cxt_maxlim = cxt_null_tot(ind_threshold)
%%
figure()
subplot(121)
histogram(cxt_null_tot,'Normalization','pdf')
hold on 
xline(cxt_maxlim,'Color','r','LineWidth',1.2)
xlabel('$|c_{xt}|$','Interpreter','latex')
ylabel('PDF')
xlim([0,0.01])
niceplot(16)
subplot(122)
h=cdfplot(cxt_null_tot);
hold on 
xline(cxt_maxlim,'Color','r','LineWidth',1.2)
h.LineWidth= 2;
xlabel('$|c_{xt}|$','Interpreter','latex')
ylabel('CDF')
xlim([0,0.01])
niceplot(16)