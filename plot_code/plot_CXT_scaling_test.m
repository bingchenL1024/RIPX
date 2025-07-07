%Bingchen Liu Mar 14, 2025
% this code is a test plot code to explore the scaling of t scale of CXT 

clear
load('/data1/bliu/data/cxt_fitanalysis_max_dxwidth.mat')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/cxt_nointerp_runinfo_G0lim0p15.mat')


get_CXT_scaling2fit

get_cxt_fitanalysis_max_dxwidth

%%
ind_nonan = ~isnan(a_tot.slp2./t_scale.slp2);
power = polyfit(log10(h.slp2(ind_nonan)),log10(a_tot.slp2(ind_nonan)./t_scale.slp2(ind_nonan)),1)
%%
figure()
subplot(121)
scatter(x_nond_tot.slp2,a_tot.slp2,50,'filled')
hold on 
scatter(x_nond_tot.slp3,a_tot.slp3,50,'filled')
hold on 
scatter(x_nond_tot.slp4,a_tot.slp4,50,'filled')
legend('slp2','slp3','slp4')
xlabel('$nondx)$','Interpreter','latex')
ylabel('${\tau (\mathrm{s})}$','Interpreter','latex')
niceplot(20)
hold off  

subplot(122)
scatter(x_nond_tot.slp2,G0_nond.slp2,50,'filled')
hold on 
scatter(x_nond_tot.slp3,G0_nond.slp3,50,'filled')
hold on 
scatter(x_nond_tot.slp4,G0_nond.slp4,50,'filled')
legend('slp2','slp3','slp4')
xlabel('$xnond$','Interpreter','latex')
ylabel('$Gnond$','Interpreter','latex')
niceplot(20)
%set(gca,'XScale','log','YScale','log')
hold off 
