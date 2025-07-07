% Bingchen Liu Oct 4, 2024
% This code plot the distribution of rsq for different fit and evaluate the
% quality of the fit 


%% c= sqrt(gh)
%load('/data1/bliu/data/cxt_alongct_x_fitpara.mat')
load('/data1/bliu/data/cxt_alongct_nointerp_fitpara')
rsq = [];
for i = 1:24
    rsq = [rsq;fitpara.slp2(i).rsq];
    rsq = [rsq;fitpara.slp3(i).rsq];
    rsq = [rsq;fitpara.slp4(i).rsq];

end 

gh_mean = mean(rsq);
gh_std = std(rsq);
nbad = length(find(rsq<0.7));

figure()
histogram(rsq,50)
xlabel('$r^{2} = 1 - \frac{\sum error}{\sum var}$','interpreter','latex')
ylabel('count')
niceplot_nobold_nomintick(18)
title(['$c = \sqrt{gh}$, ', '$\bar{r^{2}}=$', num2str(gh_mean),' $\sigma=$',num2str(gh_std)],'interpreter','latex')


%% c= maxvar
clear
load('/data1/bliu/data/cxt_alongct_x_maxvar_fitpara.mat')

rsq_maxvar = [];
for i = 1:24
    rsq_maxvar = [rsq_maxvar;fitpara.slp2(i).rsq];
    rsq_maxvar = [rsq_maxvar;fitpara.slp3(i).rsq];
    rsq_maxvar = [rsq_maxvar;fitpara.slp4(i).rsq];

end 

maxvar_mean = mean(rsq_maxvar);
maxvar_std = std(rsq_maxvar);
nbad_maxvar = length(find(rsq_maxvar<0.7));

figure()
histogram(rsq_maxvar,50)
xlabel('$r^{2} = 1 - \frac{\sum error}{\sum var}$','interpreter','latex')
ylabel('count')
niceplot_nobold_nomintick(18)
title(['c = $\sqrt{gh}$ * coef,',' $\bar{r^{2}}=$', num2str(maxvar_mean),' $\sigma=$',num2str(maxvar_std)],'interpreter','latex')

%% with cxt(dt)

clear
load('/data1/bliu/data/cxt_alongct_t_fitpara.mat')

rsq_t = [];
for i = 1:24
    rsq_t = [rsq_t;fitpara.slp2(i).rsq];
    rsq_t = [rsq_t;fitpara.slp3(i).rsq];
    rsq_t = [rsq_t;fitpara.slp4(i).rsq];

end 

t_mean = mean(rsq_t);
t_std = std(rsq_t);
nbad_t = length(find(rsq_t<0.7));

figure()
histogram(rsq_t,50)
xlabel('$r^{2} = 1 - \frac{\sum error}{\sum var}$','interpreter','latex')
ylabel('count')
niceplot_nobold_nomintick(18)
title(['$c_{xt}(dt)$',' $\bar{r^{2}}=$', num2str(t_mean),' $\sigma=$',num2str(t_std)],'interpreter','latex')
