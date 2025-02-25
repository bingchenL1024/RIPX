% Bingchen Liu Aug 15, 2024
% This code select run and cross-shore locations that has low rsq values
% and plot the fit to evaluate


clear
load('/data1/bliu/data/cxt_alongct_fitpara.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')

%% pick out low rsq run and x indices

%slope = 'slp2';
rsq_lim = 0.5;

for runind = 1:24
    A = fitpara.slp2(runind);
    xind_length = length(fitpara.slp2(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsq_lim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).ind = badfit_temp;
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).t{i}= A.t_itp_fit{badfit_temp(i)};
        
        %%%%%%%%%%%%%% NaN the original data
%             fitpara.slp2(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
%             fitpara.slp2(runind).a(badfit_temp(i)) = NaN;
%             fitpara.slp2(runind).b(badfit_temp(i)) = NaN;
%             fitpara.slp2(runind).c(badfit_temp(i)) = NaN;
        end 
    end 
    
    for i = 1:length(badfit_temp)
        cxt_data_badfit= A.cxt_data{badfit_temp(i)};
        cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
        t_badfit = A.t_itp_fit{badfit_temp(i)};
        
        figure()
        plot(t_badfit,cxt_data_badfit,'LineWidth',1.5);   
        hold on;
        plot(t_badfit,cxt_fit_badfit,'LineWidth',1.5);
        hold off;
        legend('data','fit')
        title(['run:',num2str(indbath.slp2(runind)),' run(local):',num2str(runind),' xind:',num2str(badfit_temp(i)),' rsq=',num2str(A.rsq(badfit_temp(i)))])
        clear cxt_fit_badfit cxt_data_badfit t_badfit
    end 
end 
        


for runind = 1:24
    A = fitpara.slp3(runind);
    xind_length = length(fitpara.slp3(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsq_lim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).ind = badfit_temp;
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).t{i}= A.t_itp_fit{badfit_temp(i)};
        end 
    end 
    
    for i = 1:length(badfit_temp)
        cxt_data_badfit= A.cxt_data{badfit_temp(i)};
        cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
        t_badfit = A.t_itp_fit{badfit_temp(i)};
        
        figure()
        plot(t_badfit,cxt_data_badfit,'LineWidth',1.5);   
        hold on;
        plot(t_badfit,cxt_fit_badfit,'LineWidth',1.5);
        hold off;
        legend('data','fit')
        title(['run:',num2str(indbath.slp3(runind)),' run(local):',num2str(runind),' xind:',num2str(badfit_temp(i)),' rsq=',num2str(A.rsq(badfit_temp(i)))])
        clear cxt_fit_badfit cxt_data_badfit t_badfit
    end 
end 


for runind = 1:24
    A = fitpara.slp4(runind);
    xind_length = length(fitpara.slp4(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsq_lim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).ind = badfit_temp;
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).t{i}= A.t_itp_fit{badfit_temp(i)};
        end 
    end 
    
    for i = 1:length(badfit_temp)
        cxt_data_badfit= A.cxt_data{badfit_temp(i)};
        cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
        t_badfit = A.t_itp_fit{badfit_temp(i)};
        
        figure()
        plot(t_badfit,cxt_data_badfit,'LineWidth',1.5);   
        hold on;
        plot(t_badfit,cxt_fit_badfit,'LineWidth',1.5);
        hold off;
        legend('data','fit')
        title(['run:',num2str(indbath.slp4(runind)),' run(local):',num2str(runind),' xind:',num2str(badfit_temp(i)),' rsq=',num2str(A.rsq(badfit_temp(i)))])
        clear cxt_fit_badfit cxt_data_badfit t_badfit
    end 
end 
%% re-fit to see why bad fit 



t_refit = badfit.slp3.runind_22.t{1};
cxt_data_refit = badfit.slp3.runind_22.data{1};



ftt = strcat('exp(-x/a)*cos(x/b-c)/cos(c)');
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 -20 -20];
opts.StartPoint = [0.5 0.25 0.4 ]; % beginning parameters - amp, mu, std.
opts.Upper = [20 20 20];
        
[f,gof]=fit(t_refit,cxt_data_refit,ft, opts);
cxt_fit = exp(-t_refit./f.a).*cos(t_refit/(f.b)-f.c)/(cos(f.c));

figure()
plot(t_refit,cxt_data_refit,'LineWidth',1.5);   
hold on;
plot(t_refit,cxt_fit,'LineWidth',1.5);
hold off;
legend('data','fit')
title(['rsq=',num2str(gof.rsquare)])
clear cxt_fit_badfit cxt_data_badfit t_badfit

%% get min a value
a_tot = [];
b_tot = [];
c_tot = [];
rsq_tot = [];
for i = 1:24
    a_tot= [a_tot; fitpara.slp2(i).a];
    b_tot= [b_tot; fitpara.slp2(i).b];
    c_tot= [c_tot; fitpara.slp2(i).c];
    rsq_tot= [rsq_tot; fitpara.slp2(i).rsq];

end 








