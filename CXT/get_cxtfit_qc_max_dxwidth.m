% Bingchen Liu March 3, 2025
% This code qc the cxt fit data with max dxwidth extraction method 
% it is based on the rsquare in the fit 



clear
load('/data1/bliu/data/cxt_alongct_nointerp_fitpara_max_dxwidth')
%load('/data1/bliu/data/cxt_alongct_x_maxvar_fitpara.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')

plot_badfit = 'No';
rsqlim = 0.98;
% get general stats of r^2 for all runs 

[rsq_all] = get_combine_3slp_24run(fitpara,'rsq',1);

lim_percent = (sum(rsq_all<rsqlim)/length(rsq_all))*100
%p_cdf= cdf(rsq_all);
% figure
% cdfplot(rsq_all)
% xlim([0.9,1])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% save bad fit data for re-fitting 

for runind = 1:24
    A = fitpara.slp2(runind);
    xind_length = length(fitpara.slp2(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).ind = badfit_temp;
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp2.(['runind_',num2str(indbath.slp2(runind))]).t{i}= A.t_itp{badfit_temp(i)};
        
        %%%%%%%%%%%%%% NaN the original data
%             fitpara.slp2(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
%             fitpara.slp2(runind).a(badfit_temp(i)) = NaN;
%             fitpara.slp2(runind).b(badfit_temp(i)) = NaN;
%             fitpara.slp2(runind).c(badfit_temp(i)) = NaN;
        end 
    end 
    if contains(plot_badfit,'Yes')
        for i = 1:length(badfit_temp)
            cxt_data_badfit= A.cxt_data{badfit_temp(i)};
            cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
            t_badfit = A.t_itp{badfit_temp(i)};

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
end 
        


for runind = 1:24
    A = fitpara.slp3(runind);
    xind_length = length(fitpara.slp3(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).ind = badfit_temp;
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp3.(['runind_',num2str(indbath.slp3(runind))]).t{i}= A.t_itp{badfit_temp(i)};
        
        %%%%%%%%%%%%%% NaN the original data
%             fitpara.slp3(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
%             fitpara.slp3(runind).a(badfit_temp(i)) = NaN;
%             fitpara.slp3(runind).b(badfit_temp(i)) = NaN;
%             fitpara.slp3(runind).c(badfit_temp(i)) = NaN;
        end 
    end 
    if contains(plot_badfit,'Yes')
        for i = 1:length(badfit_temp)
            cxt_data_badfit= A.cxt_data{badfit_temp(i)};
            cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
            t_badfit = A.t_itp{badfit_temp(i)};

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
end 


for runind = 1:24
    A = fitpara.slp4(runind);
    xind_length = length(fitpara.slp4(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);
    
    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).ind = badfit_temp;
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).data{i}= A.cxt_data{badfit_temp(i)};
            badfit.slp4.(['runind_',num2str(indbath.slp4(runind))]).t{i}= A.t_itp{badfit_temp(i)};
        
        %%%%%%%%%%%%%% NaN the original data
%             fitpara.slp4(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
%             fitpara.slp4(runind).a(badfit_temp(i)) = NaN;
%             fitpara.slp4(runind).b(badfit_temp(i)) = NaN;
%             fitpara.slp4(runind).c(badfit_temp(i)) = NaN;
        end 
    end
    
    if contains(plot_badfit,'Yes')
        for i = 1:length(badfit_temp)
            cxt_data_badfit= A.cxt_data{badfit_temp(i)};
            cxt_fit_badfit= A.cxt_fit{badfit_temp(i)};
            t_badfit = A.t_itp{badfit_temp(i)};

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
end 


save('/data1/bliu/data/cxt_badfit_nointerp.mat','badfit')


%% %%%%%%%%%%%%%%%%%% QC bad fit (turn on or off)

for runind = 1:24
    A = fitpara.slp2(runind);
    xind_length = length(fitpara.slp2(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);

    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        %%%%%%%%%%%%% NaN the original data
            fitpara.slp2(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
            fitpara.slp2(runind).a(badfit_temp(i)) = NaN;
            fitpara.slp2(runind).cxt_data{badfit_temp(i)}(:) = NaN;

        end 
    end 


end 

for runind = 1:24
    A = fitpara.slp3(runind);
    xind_length = length(fitpara.slp3(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);

    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        %%%%%%%%%%%%% NaN the original data
            fitpara.slp3(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
            fitpara.slp3(runind).a(badfit_temp(i)) = NaN;
            fitpara.slp3(runind).cxt_data{badfit_temp(i)}(:) = NaN;

        end 
    end 


end 

for runind = 1:24
    A = fitpara.slp4(runind);
    xind_length = length(fitpara.slp4(runind).x);
    rsq  = A.rsq;
    badfit_temp = find(rsq<rsqlim);

    if ~isempty(badfit_temp)
        for i = 1:length(badfit_temp)
        %%%%%%%%%%%%% NaN the original data
            fitpara.slp4(runind).cxt_fit{badfit_temp(i)}(:) = NaN;
            fitpara.slp4(runind).a(badfit_temp(i)) = NaN;
            fitpara.slp4(runind).cxt_data{badfit_temp(i)}(:) = NaN;

        end 
    end 


end 
        
head = 'generated in get_cxtfit_qc_max_dxwidth';
save('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced.mat','fitpara')
save('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_1mres.mat','fitpara')

