% Bingchen Liu, Aug 1, 2024
% This code fit the funcitonal form of cxt_alongct (analyzed from code get_cxt_alongct)

clear
close all 
load('/data1/bliu/data/cxt_alongct_t.mat')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/ind_of_diff_bath.mat')

domain_fit = 1:36; %0 ~3.5 s 
t_itp_fit = t_itp(domain_fit);
t_itp = 1;
g=9.81;


%ftt= strcat('exp(-x/a)*cos(x/a)');
%ftt = strcat('exp(-x/a)*cos(x/b)');
ftt = strcat('exp(-x/a)*cos(x/b-c)/cos(c)');
%ftt = strcat('exp(-x/a +c*cos(x/b)-c)*cos(x/b)');
%ftt = strcat('exp(-x/a +c*cos(x/b+d)-c*cos(d))*cos(x/b+d)/(cos(d))');

ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 -20 -20];
opts.StartPoint = [0.5 0.25 0.4 ]; % beginning parameters - amp, mu, std.
opts.Upper = [20 20 20];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp2


for ind_slp2 = 1:24
    runnum = indbath.slp2(ind_slp2);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_t_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;
    
%     fitpara.slp2{ind_slp2,1}.a = zeros(length(1:5:length(ind_good)),1); %get fit every 5m
%     fitpara.slp2{ind_slp2,1}.b = zeros(length(1:5:length(ind_good)),1); %get fit every 5m
%     fitpara.slp2{ind_slp2,1}.c = zeros(length(1:5:length(ind_good)),1); %get fit every 5m
%     fitpara.slp2{ind_slp2,1}.x = zeros(length(1:5:length(ind_good)),1); %get fit every 5m
%     fitpara.slp2{ind_slp2,1}.x_nond = zeros(length(1:5:length(ind_good)),1); %get fit every 5m

    
    %cxt_alongct_1run = cxt_alongct_mean{runnum};
    %cxt_alongct_1run = cxt_alongct_1run(:,70);
    
    
    
    j = 0;
    for i = 1:5:dim_cxt(2) %crossshore scale ----------fit every 5 meter
        j = j+1;
        cxt_data = cxt_alongct_1run(:,i);
        cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
        [f,gof]=fit(t_itp_fit(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp_fit(include)./f.a).*cos(t_itp_fit(include)/(f.b)-f.c)/(cos(f.c));
%%%%%%%%%%%%%% Quick QC
% figure(1);plot(t_itp_fit(include)',cxt_data_fit(include)); hold on;plot(t_itp_fit(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)
        fitpara.slp2{ind_slp2,1}.x(j,1) = x(i);        
        fitpara.slp2{ind_slp2,1}.x_nond(j,1) = x_nond(i); 
        fitpara.slp2{ind_slp2,1}.xb(j,1) = xb;        
        fitpara.slp2{ind_slp2,1}.t_scale(j,1) = t_scale(i);         
        fitpara.slp2{ind_slp2,1}.a(j,1) = f.a;
        fitpara.slp2{ind_slp2,1}.b(j,1) = f.b;
        fitpara.slp2{ind_slp2,1}.c(j,1) = f.c; 
        fitpara.slp2{ind_slp2,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp2{ind_slp2,1}.t_itp_fit{j,1} = t_itp_fit(include)';        
        fitpara.slp2{ind_slp2,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp2{ind_slp2,1}.cxt_fit{j,1} = cxt_fit';        
   
        clear f cxt_data_fit cxt_fit
        
        
%%%%%%%%%%%%%%%%%%%%%%%\\\\\\\\\\\\\\\\\\\\\\ for QC
% 
%         %cxt_fit = exp(-t_itp(include)./f.a).*cos(t_itp(include)./f.a);
%         %cxt_fit = exp(-t_itp./f.a).*cos(t_itp/(f.b));
%         cxt_fit = exp(-t_itp_fit(include)./f.a).*cos(t_itp_fit(include)/(f.b)-f.c)/(cos(f.c));
%         %cxt_fit = exp(-t_itp(include)/f.a +f.c*cos(t_itp(include)/f.b)-f.c).*cos(t_itp(include)/(f.b));
%         %cxt_fit = exp(-t_itp(include)/f.a +f.c*cos(t_itp(include)/f.b+f.d)-f.c*(cos(f.d))).*cos(t_itp(include)/(f.b)+f.d)/(cos(f.d));
% 
%         misfit = cxt_data_fit(include)-cxt_fit';
%        
%         
%         figure()
%         plot(t_itp_fit(include),cxt_data_fit(include),'k','LineWidth',2)
%         hold on 
%         plot(t_itp_fit(include),cxt_data_fit(include),'k','LineWidth',2)
%         hold on 
%         plot(t_itp_fit(include),cxt_fit,'b','LineWidth',2)
%         hold on 
%         plot(t_itp_fit(include),misfit,'r','LineWidth',2)        
%         hold off 
%         legend('data(whole domain)','data(fit domain)','fit','misfit')
%         xlim([0,4])
%         niceplot_nobold_nomintick(18);
%         xlabel('Time Lag (s)--along x=ct')
%         ylabel('Cross Correlation')
%         

%%%%%%%%%%%%%%%%%%%%\\\\\\\\\\\\\\\\\\\\\\\ fit only the exponential (for QC)
    %%%%%%%%%%%%%%%%%%%%%%5 fit only the exponential 
%       ftt_exp = strcat('exp(-(x/a)+c*cos(x/b)-c)');
%       ft_exp = fittype( sprintf('%s',ftt_exp));
%        test_x= [0, 0.5 1.2 1.7 2.2 2.7 3.2];
%        test_y= [1 0.158 0.1289 0.041 0.0498 0.01 0.028];
%        [f_exp,gof_exp]=fit(test_x',test_y',ft_exp);
%        test_y_fit =exp(-(test_x/f.a)+f.c*cos(test_x/f.b)-f.c);
% 
% 
%         
%         
%        
%         
%         figure()
%         plot(test_x,test_y,'k','LineWidth',2)
%         hold on 
%         plot(test_x,test_y_fit,'b','LineWidth',2)
%         legend('data','fit','misfit')
%         hold off 
        
    end 


ind_slp2
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp3

for ind_slp3 = 1:24
    runnum = indbath.slp3(ind_slp3);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_t_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;


    
    j = 0;
    for i = 1:5:dim_cxt(2) %fit every 5 meter
        j = j+1;
        cxt_data = cxt_alongct_1run(:,i);
        cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
%         figure()
%         plot(t_itp_fit(include)',cxt_fit)
%         hold on 
%         plot(t_itp_fit(include)',cxt_data_fit(include))
%         hold off 
        [f,gof]=fit(t_itp_fit(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp_fit(include)./f.a).*cos(t_itp_fit(include)/(f.b)-f.c)/(cos(f.c));
%figure(1);plot(t_itp_fit(include)',cxt_data_fit(include)); hold on;plot(t_itp_fit(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)

        
        fitpara.slp3{ind_slp3,1}.x(j,1) = x(i);        
        fitpara.slp3{ind_slp3,1}.x_nond(j,1) = x_nond(i);  
        fitpara.slp3{ind_slp3,1}.xb(j,1) = xb;        
        fitpara.slp3{ind_slp3,1}.t_scale(j,1) = t_scale(i);         
        fitpara.slp3{ind_slp3,1}.a(j,1) = f.a;
        fitpara.slp3{ind_slp3,1}.b(j,1) = f.b;
        fitpara.slp3{ind_slp3,1}.c(j,1) = f.c; 
        fitpara.slp3{ind_slp3,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp3{ind_slp3,1}.t_itp_fit{j,1} = t_itp_fit(include)';        
        fitpara.slp3{ind_slp3,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp3{ind_slp3,1}.cxt_fit{j,1} = cxt_fit';        
   
    clear f cxt_data_fit cxt_fit
        
    end 
ind_slp3
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp4
for ind_slp4 = 1:24
    runnum = indbath.slp4(ind_slp4);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_t_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;
    
    j = 0;
    for i = 1:5:dim_cxt(2) %fit every 5 meter
        j = j+1;
        cxt_data = cxt_alongct_1run(:,i);
        cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
%         figure()
%         plot(t_itp(include)',cxt_data(include))
        [f,gof]=fit(t_itp_fit(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp_fit(include)./f.a).*cos(t_itp_fit(include)/(f.b)-f.c)/(cos(f.c));
% figure(1);plot(t_itp_fit(include)',cxt_data_fit(include)); hold on;plot(t_itp_fit(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)

        fitpara.slp4{ind_slp4,1}.x(j,1) = x(i);        
        fitpara.slp4{ind_slp4,1}.x_nond(j,1) = x_nond(i);  
        fitpara.slp4{ind_slp4,1}.xb(j,1) = xb;        
        fitpara.slp4{ind_slp4,1}.t_scale(j,1) = t_scale(i);         
        fitpara.slp4{ind_slp4,1}.a(j,1) = f.a;
        fitpara.slp4{ind_slp4,1}.b(j,1) = f.b;
        fitpara.slp4{ind_slp4,1}.c(j,1) = f.c; 
        fitpara.slp4{ind_slp4,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp4{ind_slp4,1}.t_itp_fit{j,1} = t_itp_fit(include)';        
        fitpara.slp4{ind_slp4,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp4{ind_slp4,1}.cxt_fit{j,1} = cxt_fit';        
   
    clear f cxt_data_fit cxt_fit
        
    end 
ind_slp4
end



%% convert it to be more compact

fitpara.slp2  = cell2mat(fitpara.slp2);
fitpara.slp3  = cell2mat(fitpara.slp3);
fitpara.slp4  = cell2mat(fitpara.slp4);

readme = 'created in get_cxt_fit';

save('/data1/bliu/data/cxt_alongct_t_fitpara','fitpara','readme')