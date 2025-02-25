% Bingchen Liu Nov 13, 2024
% This code fit cxt along ct without interpolation 



clear
close all 
load('/data1/bliu/data/cxt_alongct_nointerp.mat')
%load('/data1/bliu/data/cxt_alongct_x_maxvar.mat')
load('/data1/bliu/data/SS_raw.mat') 
load('/data1/bliu/data/ind_of_diff_bath.mat')
load('/data1/bliu/data/cxt_ind_good.mat')

%%

%domain_fit = 1:length(x_itp); %0 ~3.5 s 
%t_itp = x_itp(domain_fit);

xskip = 5;

g=9.81;


%ftt= strcat('exp(-x/a)*cos(x/a)');
%ftt = strcat('exp(-x/a)*cos(x/b)');
ftt = strcat('exp(-x/a)');
%ftt = strcat('exp(-x/a +c*cos(x/b)-c)*cos(x/b)');
%ftt = strcat('exp(-x/a +c*cos(x/b+d)-c*cos(d))*cos(x/b+d)/(cos(d))');

ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 0 -pi];
opts.StartPoint = [0.5]; % beginning parameters - amp, mu, std.
%opts.Upper = [20 20 pi]; %phase shift should be within 2pi


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp2


for ind_slp2 = 1:24
    runnum = indbath.slp2(ind_slp2);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    if_interp = if_interp_All{runnum};%ifinterp
    
    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification
    
    j = 0;
    for i = 1:xskip:dim_cxt(2) %crossshore scale ----------fit every 5 meter
        j = j+1;
        cxt_data_fit = cxt_alongct_1run(:,i);
        %cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
        [f,gof]=fit(t_itp(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp(include)./f.a);
%%%%%%%%%%%%%% Quick QC
% figure(1);plot(t_itp(include)',cxt_data_fit(include)); hold on;plot(t_itp(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)
        fitpara.slp2{ind_slp2,1}.x(j,1) = x(i);        
        fitpara.slp2{ind_slp2,1}.x_nond(j,1) = x_nond(i); 
        fitpara.slp2{ind_slp2,1}.G0_nond(j,1) = G0_nond(i); 
        fitpara.slp2{ind_slp2,1}.xb(j,1) = xb;  
        fitpara.slp2{ind_slp2,1}.Tp(j,1) = Tp;        
        fitpara.slp2{ind_slp2,1}.t_scale(j,1) = t_scale(i);    
        fitpara.slp2{ind_slp2,1}.h(j,1) = h(i);  
        fitpara.slp2{ind_slp2,1}.kw(j,1) = kw(i);  

        %fitpara.slp2{ind_slp2,1}.Hs(j,1) = Hs(i);          
        fitpara.slp2{ind_slp2,1}.a(j,1) = f.a;
        fitpara.slp2{ind_slp2,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp2{ind_slp2,1}.t_itp{j,1} = t_itp(include)';        
        fitpara.slp2{ind_slp2,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp2{ind_slp2,1}.cxt_fit{j,1} = cxt_fit';        
        fitpara.slp2{ind_slp2,1}.runnum{j,1} = runnum;           
        clear f cxt_data_fit cxt_fit
        
        

    end 
fitpara.slp2{ind_slp2,1}.Hs_interp = interp1(SS.X2,SS.Hs,fitpara.slp2{ind_slp2,1}.x); 


%ind_slp2
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp3

for ind_slp3 = 1:24
    runnum = indbath.slp3(ind_slp3);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    if_interp = if_interp_All{runnum};%ifinterp

    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification


    
    
    j = 0;
    for i = 1:5:dim_cxt(2) %fit every 5 meter
        j = j+1;
        cxt_data_fit = cxt_alongct_1run(:,i);
        %cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
%         figure()
%         plot(t_itp(include)',cxt_fit)
%         hold on 
%         plot(t_itp(include)',cxt_data_fit(include))
%         hold off 
        [f,gof]=fit(t_itp(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp(include)./f.a);
%figure(1);plot(t_itp(include)',cxt_data_fit(include)); hold on;plot(t_itp(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)

        
        fitpara.slp3{ind_slp3,1}.x(j,1) = x(i);        
        fitpara.slp3{ind_slp3,1}.x_nond(j,1) = x_nond(i);  
        fitpara.slp3{ind_slp3,1}.G0_nond(j,1) = G0_nond(i); 
        fitpara.slp3{ind_slp3,1}.xb(j,1) = xb;  
        fitpara.slp3{ind_slp3,1}.Tp(j,1) = Tp;        
        fitpara.slp3{ind_slp3,1}.t_scale(j,1) = t_scale(i);   
        fitpara.slp3{ind_slp3,1}.h(j,1) = h(i);   
        fitpara.slp3{ind_slp3,1}.kw(j,1) = kw(i);  
        fitpara.slp3{ind_slp3,1}.a(j,1) = f.a;
        fitpara.slp3{ind_slp3,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp3{ind_slp3,1}.t_itp{j,1} = t_itp(include)';        
        fitpara.slp3{ind_slp3,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp3{ind_slp3,1}.cxt_fit{j,1} = cxt_fit';        
        fitpara.slp3{ind_slp3,1}.runnum{j,1} = runnum;        

    clear f cxt_data_fit cxt_fit
        
    end 
fitpara.slp3{ind_slp3,1}.Hs_interp = interp1(SS.X2,SS.Hs,fitpara.slp3{ind_slp3,1}.x); 

%ind_slp3
end 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% slp4
for ind_slp4 = 1:24
    runnum = indbath.slp4(ind_slp4);
    ind_good = ind_good_All{runnum};
    cxt_alongct_1run  = cxt_alongct_ALL{runnum};
    dim_cxt= size(cxt_alongct_1run);
    if_interp = if_interp_All{runnum};%ifinterp


    SS = S(runnum);
    x = SS.X(ind_good);
    x_nond = x/SS.xb;
    h = SS.h(ind_good);
    t_scale = sqrt(h/g);
    xb = SS.xb;
    Tp = SS.Tp_T;
    kw = get_wavenum(2*pi/(Tp),h);
    G0_nond = abs(SS.curlF_std(ind_good))./abs(max(SS.curlF_std(ind_good)));%G0 modification

    
    
    j = 0;
    for i = 1:5:dim_cxt(2) %fit every 5 meter
        j = j+1;
        cxt_data_fit = cxt_alongct_1run(:,i);
        %cxt_data_fit = cxt_alongct_1run(domain_fit,i);
        include = ~isnan(cxt_data_fit);
%         figure()
%         plot(t_itp(include)',cxt_data(include))
        [f,gof]=fit(t_itp(include)',cxt_data_fit(include),ft, opts);
        cxt_fit = exp(-t_itp(include)./f.a);
% figure(1);plot(t_itp(include)',cxt_data_fit(include)); hold on;plot(t_itp(include)',cxt_fit);hold off;legend('data','fit')
% pause(0.2)

        fitpara.slp4{ind_slp4,1}.x(j,1) = x(i);        
        fitpara.slp4{ind_slp4,1}.x_nond(j,1) = x_nond(i);
        fitpara.slp4{ind_slp4,1}.G0_nond(j,1) = G0_nond(i); 
        fitpara.slp4{ind_slp4,1}.xb(j,1) = xb; 
        fitpara.slp4{ind_slp4,1}.Tp(j,1) = Tp;                
        fitpara.slp4{ind_slp4,1}.t_scale(j,1) = t_scale(i);
        fitpara.slp4{ind_slp4,1}.h(j,1) = h(i);
        fitpara.slp4{ind_slp4,1}.kw(j,1) = kw(i);  
        fitpara.slp4{ind_slp4,1}.a(j,1) = f.a; 
        fitpara.slp4{ind_slp4,1}.rsq(j,1) = gof.rsquare;        
        fitpara.slp4{ind_slp4,1}.t_itp{j,1} = t_itp(include)';        
        fitpara.slp4{ind_slp4,1}.cxt_data{j,1} = cxt_data_fit(include);        
        fitpara.slp4{ind_slp4,1}.cxt_fit{j,1} = cxt_fit';     
        fitpara.slp4{ind_slp4,1}.runnum{j,1} = runnum;        

   
    clear f cxt_data_fit cxt_fit
        
    end 
fitpara.slp4{ind_slp4,1}.Hs_interp = interp1(SS.X2,SS.Hs,fitpara.slp4{ind_slp4,1}.x); 

%ind_slp4
end



%% convert it to be more compact

fitpara.slp2  = cell2mat(fitpara.slp2);
fitpara.slp3  = cell2mat(fitpara.slp3);
fitpara.slp4  = cell2mat(fitpara.slp4);

readme = 'created in get_cxt_fit_nointerp (fit cxt(dx))';
save('/data1/bliu/data/cxt_alongct_nointerp_fitpara','fitpara','readme')
