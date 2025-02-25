% Bingchen Liu Aug 15, 2024
% This code analysis a b c fit from fit of cxt 

clear
load('/data1/bliu/data/cxt_alongct_t_fitpara.mat')
g=9.81;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp2
a_tot.slp2 = [];
b_tot.slp2 = [];
c_tot.slp2 = [];
rsq_tot.slp2 = [];
x_tot.slp2 = [];
x_br_tot.slp2 = [];
x_nond_tot.slp2 = [];
t_scale.slp2=[];
%h.slp2 = [];
for i = 1:24
    a_tot.slp2= [a_tot.slp2; fitpara.slp2(i).a];
    b_tot.slp2= [b_tot.slp2; fitpara.slp2(i).b];
    c_tot.slp2= [c_tot.slp2; fitpara.slp2(i).c];
    x_tot.slp2 = [x_tot.slp2;fitpara.slp2(i).x];
    x_br_tot.slp2 = [x_br_tot.slp2;fitpara.slp2(i).xb];
    x_nond_tot.slp2 = [x_nond_tot.slp2; fitpara.slp2(i).x_nond];
    rsq_tot.slp2= [rsq_tot.slp2; fitpara.slp2(i).rsq];
    t_scale.slp2=[t_scale.slp2;fitpara.slp2(i).t_scale];
    %h.slp2=[h.slp2;fitpara.slp2(i).h];
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
%h.slp3 = [];

for i = 1:24
    a_tot.slp3= [a_tot.slp3; fitpara.slp3(i).a];
    b_tot.slp3= [b_tot.slp3; fitpara.slp3(i).b];
    c_tot.slp3= [c_tot.slp3; fitpara.slp3(i).c];
    x_tot.slp3 = [x_tot.slp3;fitpara.slp3(i).x];
    x_br_tot.slp3 = [x_br_tot.slp3;fitpara.slp3(i).xb];
    x_nond_tot.slp3 = [x_nond_tot.slp3; fitpara.slp3(i).x_nond];
    rsq_tot.slp3= [rsq_tot.slp3; fitpara.slp3(i).rsq];
    t_scale.slp3=[t_scale.slp3;fitpara.slp3(i).t_scale];
    %h.slp3=[h.slp3;fitpara.slp3(i).h];

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
%h.slp4 = [];

for i = 1:24
    a_tot.slp4= [a_tot.slp4; fitpara.slp4(i).a];
    b_tot.slp4= [b_tot.slp4; fitpara.slp4(i).b];
    c_tot.slp4= [c_tot.slp4; fitpara.slp4(i).c];
    x_tot.slp4 = [x_tot.slp4;fitpara.slp4(i).x];
    x_br_tot.slp4 = [x_br_tot.slp4;fitpara.slp4(i).xb];
    x_nond_tot.slp4 = [x_nond_tot.slp4; fitpara.slp4(i).x_nond];
    rsq_tot.slp4= [rsq_tot.slp4; fitpara.slp4(i).rsq];
    t_scale.slp4=[t_scale.slp4;fitpara.slp4(i).t_scale];
    %h.slp4=[h.slp4;fitpara.slp4(i).h];

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

%% a and b vs x 

figure()
scatter(x_nond_tot.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(x_nond_tot.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(x_nond_tot.slp4,a_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Dimensionless surfzone location (x/xb)')
ylabel('Fit parameter a (decay rate)')
legend('slp2','slp3','slp4')
xlim([-1,0])
hold off 

figure()
scatter(x_nond_tot.slp2,b_tot.slp2,45,'filled')
hold on 
scatter(x_nond_tot.slp3,b_tot.slp3,45,'filled')
hold on 
scatter(x_nond_tot.slp4,b_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Dimensionless surfzone location (x/xb)')
ylabel('Fit parameter b (oscillation freq)')
legend('slp2','slp3','slp4')
xlim([-1,0])
ylim([-inf 0.5])
hold off 

figure()
scatter(x_nond_tot.slp2,c_tot.slp2,45,'filled')
hold on 
scatter(x_nond_tot.slp3,c_tot.slp3,45,'filled')
hold on 
scatter(x_nond_tot.slp4,c_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Dimensionless surfzone location (x/xb)')
ylabel('Fit parameter c (phase shift)')
legend('slp2','slp3','slp4')
hold off 
%% a vs t scale

figure()
scatter(t_scale.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(t_scale.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(t_scale.slp4,a_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Time scale ($\sqrt{h/g}$)','interpreter','latex')
ylabel('Fit parameter a (decay rate)')
legend('slp2','slp3','slp4')
hold off 

figure()
scatter(-t_scale.slp2.*x_nond_tot.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(-t_scale.slp3.*x_nond_tot.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(-t_scale.slp4.*x_nond_tot.slp4,a_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Time scale ($\sqrt{h/g} \, \frac{x}{x_{b}}$)','interpreter','latex')
ylabel('Fit parameter a (decay rate)')
legend('slp2','slp3','slp4')
hold off 

figure()
scatter(t_sincebr.slp2,a_tot.slp2,45,'filled')
hold on 
scatter(t_sincebr.slp3,a_tot.slp3,45,'filled')
hold on 
scatter(t_sincebr.slp4,a_tot.slp4,45,'filled')
niceplot_nobold_nomintick(18)
xlabel('Time since breaking','interpreter','latex')
ylabel('Fit parameter a (decay rate)')
legend('slp2','slp3','slp4')
hold off 
%% relation between a and b 

include = find(a_tot.slp2>0.1 & a_tot.slp2<1&b_tot.slp2<0.5&b_tot.slp2>0);

%ftt = strcat('real(acos(a*exp(1/(x))))') ;
ftt = strcat('real(b*acos(a*exp(1/(x)))+c)') ;
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-20 -20 -20];
opts.StartPoint = [-0.1 0 0]; % beginning parameters - amp, mu, std.
opts.Upper = [0 20 20];
[f,gof]=fit(a_tot.slp2(include),b_tot.slp2(include),ft, opts);
a_fit.slp2= 0:0.1:1.5;
b_fit.slp2 = real(f.b.*acos(f.a.*exp(1./(a_fit.slp2)))+f.c);
%b_fit = real(acos(f.a.*exp(1./a_fit)));


ftt1 = strcat('1/(c*(x+a))+b');
ft1 = fittype( sprintf('%s',ftt1));
opts1 = fitoptions( ft1 );
opts1.Display = 'Off';
opts1.Lower = [-inf -20 -inf];
opts1.StartPoint = [-0.1 0 1]; % beginning parameters - amp, mu, std.
opts1.Upper = [-20 20 inf];
[f1,gof1]=fit(a_tot.slp2(include),b_tot.slp2(include),ft1, opts1);
a_fit.slp2= 0:0.1:1.5;
b_fit1.slp2 = (1./(f1.c.*(a_fit.slp2+f1.a))) + f1.b;

figure()
scatter(a_tot.slp2(include),b_tot.slp2(include),50,'filled','b')
hold on 
plot(a_fit.slp2,b_fit.slp2,'LineWidth',2)
%plot(a_fit.slp2,b_fit1.slp2,'LineWidth',2)
xlabel('a')
ylabel('b')
niceplot_nobold_nomintick(18)
legend('data','fit')
title('Using fit = real(b*acos(a*exp(1/(x)))+c)')
hold off 
%ylim([0 0.5])


