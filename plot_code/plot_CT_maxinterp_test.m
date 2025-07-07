% Bingchen Liu May 13, 2025
% This code plot c_T vs dt that uses cxt interp on c_{fit} velocity that is
% the fit of locations of max cxt

clear
load('/data1/bliu/data/CXT_all')
load('/data1/bliu/data/ind_of_diff_bath.mat')


ind_slp234 = [indbath.slp2;indbath.slp3;indbath.slp4];

for i = 1:72
    ct_interp{i} = cxt_alongct_itp_ALL{ind_slp234(i)};
    ct_nointerp{i}= cxt_alongct_ALL{ind_slp234(i)};
end 

figure()
subplot(121)
for runind = 1:72
    plot(0:5,ct_interp{runind},'LineWidth',1)
    hold on 
end 
hold off 
xlabel('$\Delta t$','Interpreter','latex','FontSize',16)
ylabel('$c_{XT}$','Interpreter','latex','FontSize',16)
title('CXT extracted using interpolation from c_{fit}')
ylim([-0.2,1])

subplot(122)
for runind = 1:72
    plot(0:5,ct_nointerp{runind},'LineWidth',1)
    hold on 
end 
hold off 
xlabel('$\Delta t$','Interpreter','latex','FontSize',16)
ylabel('$c_{XT}$','Interpreter','latex','FontSize',16)
title('CXT extracted using max cxt method')