% Bingchen Liu Oct 7, 2024
% This code analyze the max var velocity c = sqrt(gh)*c_coef vs sqrt(gh)

load('/data1/bliu/data/cxt_alongct_x_maxvar')

c_maxvar= [];
for i = 1:120
    c_maxvar = [c_maxvar,c_coef_maxvar_All{i}];
    l(i)= length(c_coef_maxvar_All{i});
end 

figure(1)
histogram(c_maxvar,15)
xlabel('$c_{coef}$','interpreter','latex','FontSize',18)
ylabel('count','FontSize',18)
title('Histogram of c_{coef}')
niceplot(18)
