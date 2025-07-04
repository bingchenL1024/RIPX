% Bingchen Liu Mar 3, 2025
% This code calculate the nond dx width and to for plots of CXT vs dx width
% and nond CXT vs nond dx width

%dx_nond1 = dx/(c \tao(deco t scale))
%dx_nond2 = dx/(c dt(time lag))

clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth')
load('/data1/bliu/data/runnum_72run.mat')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced.mat')
load('/data1/bliu/data/ind_of_diff_bath.mat')
load('/data1/bliu/data/cxt_alongct_max_dxwidth_fitpara_qced_5loc')

%%
dx_width = -4:1:4;
g= 9.81;
alpha = 1/3;

for ind = 1:24
    runnum = indbath.slp2(ind);
    a_decotscale = fitpara.slp2(ind).a;
    for xloc = 1:length(a_decotscale)
        %cxt_dxwidth.slp2{ind,1}.dx_nond1{1,xloc} = dx_width/((fitpara.slp2(ind).h(xloc)))^alpha;
        cxt_dxwidth.slp2{ind,1}.dx_nond1{1,xloc} = dx_width/(c_fit{runnum}(xloc));
        cxt_dxwidth.slp2{ind,1}.dx_nond2{1,xloc} = dx_width/(c_fit{runnum}(xloc)*fitpara.slp2(ind).a(xloc));
        for tlag = 1:3
            cxt_dxwidth.slp2{ind,1}.cxt_nond{tlag,xloc} = cxt_pm2dx_nond_ALL{runnum}{tlag,xloc};
            cxt_dxwidth.slp2{ind,1}.cxt{tlag,xloc} = cxt_pm2dx_ALL{runnum}{tlag,xloc};

        end 

    end 
    cxt_dxwidth.slp2{ind,1}.x_nond = x_nond_All{runnum};
end 


for ind = 1:24
    runnum = indbath.slp3(ind);
    a_decotscale = fitpara.slp3(ind).a;
    for xloc = 1:length(a_decotscale)
        %cxt_dxwidth.slp3{ind,1}.dx_nond1{1,xloc} = dx_width/((fitpara.slp2(ind).h(xloc)))^alpha;
        cxt_dxwidth.slp3{ind,1}.dx_nond1{1,xloc} = dx_width/(c_fit{runnum}(xloc));
        cxt_dxwidth.slp3{ind,1}.dx_nond2{1,xloc} = dx_width/(c_fit{runnum}(xloc)*fitpara.slp3(ind).a(xloc));
        for tlag = 1:3
            cxt_dxwidth.slp3{ind,1}.cxt_nond{tlag,xloc} = cxt_pm2dx_nond_ALL{runnum}{tlag,xloc};
            cxt_dxwidth.slp3{ind,1}.cxt{tlag,xloc} = cxt_pm2dx_ALL{runnum}{tlag,xloc};
        end 
    end 
    cxt_dxwidth.slp3{ind,1}.x_nond = x_nond_All{runnum};

end 

for ind = 1:24
    runnum = indbath.slp4(ind);
    a_decotscale = fitpara.slp4(ind).a;
    for xloc = 1:length(a_decotscale)
        %cxt_dxwidth.slp4{ind,1}.dx_nond1{1,xloc} = dx_width/((fitpara.slp2(ind).h(xloc)))^alpha;
        cxt_dxwidth.slp4{ind,1}.dx_nond1{1,xloc} = dx_width/(c_fit{runnum}(xloc));
        cxt_dxwidth.slp4{ind,1}.dx_nond2{1,xloc} = dx_width/(c_fit{runnum}(xloc)*fitpara.slp4(ind).a(xloc));
        for tlag = 1:3
            cxt_dxwidth.slp4{ind,1}.cxt_nond{tlag,xloc} = cxt_pm2dx_nond_ALL{runnum}{tlag,xloc};
            cxt_dxwidth.slp4{ind,1}.cxt{tlag,xloc} = cxt_pm2dx_ALL{runnum}{tlag,xloc};
        end 
    end 
    cxt_dxwidth.slp4{ind,1}.x_nond = x_nond_All{runnum};

end 


cxt_dxwidth.slp2  = cell2mat(cxt_dxwidth.slp2);
cxt_dxwidth.slp3  = cell2mat(cxt_dxwidth.slp3);
cxt_dxwidth.slp4  = cell2mat(cxt_dxwidth.slp4);


%% get 5 loc 


loc5 = linspace(-0.75,-0.25,5);

for i = 1:24
    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(cxt_dxwidth.slp2(i).x_nond - loc5(j))); %ind out of 87*1
    
        cxt_dxwidth_5loc.slp2(i).x_nond(j) = cxt_dxwidth.slp2(i).x_nond(ind_5loc_temp(j));
        cxt_dxwidth_5loc.slp2(i).dx_nond1(j,:) = cxt_dxwidth.slp2(i).dx_nond1{ind_5loc_temp(j)}; %note: diff row rep diff xnond
        cxt_dxwidth_5loc.slp2(i).dx_nond2(j,:) = cxt_dxwidth.slp2(i).dx_nond2{ind_5loc_temp(j)};
        for tlag = 1:3
            cxt_dxwidth_5loc.slp2(i).cxt_nond{tlag}(j,:) = cxt_dxwidth.slp2(i).cxt_nond{tlag,ind_5loc_temp(j)};
        end 

    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(cxt_dxwidth.slp3(i).x_nond - loc5(j)));
    
        cxt_dxwidth_5loc.slp3(i).x_nond(j) = cxt_dxwidth.slp3(i).x_nond(ind_5loc_temp(j));
        cxt_dxwidth_5loc.slp3(i).dx_nond1(j,:) = cxt_dxwidth.slp3(i).dx_nond1{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp3(i).dx_nond2(j,:) = cxt_dxwidth.slp3(i).dx_nond2{ind_5loc_temp(j)};
        for tlag = 1:3
            cxt_dxwidth_5loc.slp3(i).cxt_nond{tlag}(j,:) = cxt_dxwidth.slp3(i).cxt_nond{tlag,ind_5loc_temp(j)};
        end 
    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(cxt_dxwidth.slp4(i).x_nond - loc5(j)));
    
        cxt_dxwidth_5loc.slp4(i).x_nond(j) = cxt_dxwidth.slp4(i).x_nond(ind_5loc_temp(j));
        cxt_dxwidth_5loc.slp4(i).dx_nond1(j,:) = cxt_dxwidth.slp4(i).dx_nond1{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp4(i).dx_nond2(j,:) = cxt_dxwidth.slp4(i).dx_nond2{ind_5loc_temp(j)};
        for tlag = 1:3
            cxt_dxwidth_5loc.slp4(i).cxt_nond{tlag}(j,:) = cxt_dxwidth.slp4(i).cxt_nond{tlag,ind_5loc_temp(j)};
        end 
    end 



end 

%% take avg/binmean of dx structure at 5 loc 


% dxlag =5:9;
% cxt_dxwidth_all = [];
% for ind = 1:24
%         cxt_temp = [cxt_dxwidth_5loc.slp2(ind).cxt_nond{1,1}(:,dxlag);cxt_dxwidth_5loc.slp3(ind).cxt_nond{1,1}(:,dxlag);cxt_dxwidth_5loc.slp4(ind).cxt_nond{1,1}(:,dxlag)];
%         cxt_dxwidth_all = [cxt_dxwidth_all;cxt_temp];
% 
% end 
dxlag =5:9;
cxt_dxwidth_5loc_all = [];
dx_nond1_5loc_all = [];
dx_nond2_5loc_all = [];

for ind = 1:24
    cxt_temp = [cxt_dxwidth_5loc.slp2(ind).cxt_nond{1,1}(:,dxlag);cxt_dxwidth_5loc.slp3(ind).cxt_nond{1,1}(:,dxlag);cxt_dxwidth_5loc.slp4(ind).cxt_nond{1,1}(:,dxlag)];
    cxt_dxwidth_5loc_all = [cxt_dxwidth_5loc_all;cxt_temp(:)];

    dx_temp1 = [cxt_dxwidth_5loc.slp2(ind).dx_nond1(:,dxlag);cxt_dxwidth_5loc.slp3(ind).dx_nond1(:,dxlag);cxt_dxwidth_5loc.slp4(ind).dx_nond1(:,dxlag)];
    dx_nond1_5loc_all = [dx_nond1_5loc_all;dx_temp1(:)];    
    
    dx_temp2 = [cxt_dxwidth_5loc.slp2(ind).dx_nond2(:,dxlag);cxt_dxwidth_5loc.slp3(ind).dx_nond2(:,dxlag);cxt_dxwidth_5loc.slp4(ind).dx_nond2(:,dxlag)];
    dx_nond2_5loc_all = [dx_nond2_5loc_all;dx_temp2(:)];     
end 


binnum = 10;
dx1lim = 1.2;
dx2lim = 2.0;

binmean1 = cxt_binmean(cxt_dxwidth_5loc_all,dx_nond1_5loc_all,binnum,dx1lim);
binmean2 = cxt_binmean(cxt_dxwidth_5loc_all,dx_nond2_5loc_all,binnum,dx2lim);

binmean2.cxt_binmean = [1,binmean2.cxt_binmean]; %add (0,1) point 
binmean2.error_binmean = [0,binmean2.error_binmean]; %add (0,1) point 
binmean2.dx_bincenter = [0,binmean2.dx_bincenter]; %add (0,1) point 

%% fit
close all

% % =================> exp fit 
% include = find(~isnan(dx_nond2_5loc_all)&~isnan(cxt_dxwidth_5loc_all)&cxt_dxwidth_5loc_all~=1);
% ftt = strcat('exp(-x/a)*cos(x/b+c)/cos(c)');
% ft = fittype( sprintf('%s',ftt));
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.StartPoint = [0.29, 0.56, 1.27]; % beginning parameters - amp, mu, std.
% %opts.Lower = [0.001,-inf,-inf];
% [fit.f,fit.gof]=fit(dx_nond2_5loc_all(include),cxt_dxwidth_5loc_all(include),ft, opts);
% 
% %[pfit_p,pfit_s] = polyfitB(dx_nond2_5loc_all(include),cxt_dxwidth_5loc_all(include),5,1);
% fit.dx_fit = 0:0.01:3;
% %cxt_fit = sinc(f.a*dx_fit);
% fit.cxt_fit = exp(-fit.dx_fit/fit.f.a).*cos(fit.dx_fit/fit.f.b+fit.f.c)/cos(fit.f.c);
% %cxt_pfit = polyval(pfit_p,dx_fit);


% ===================> R^x fit 
include = find(~isnan(dx_nond2_5loc_all)&~isnan(cxt_dxwidth_5loc_all)&cxt_dxwidth_5loc_all~=1);
ftt= strcat('a^x*cos(b*x-c)/cos(c)');
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
opts.StartPoint = [0.03, 1.785, -1.27]; % beginning parameters - amp, mu, std.
opts.Lower = [0.001,-inf,-inf];
[fit.f,fit.gof]=fit(dx_nond2_5loc_all(include),cxt_dxwidth_5loc_all(include),ft, opts);

fit.dx_fit = 0:0.01:3;
fit.cxt_fit = exp(-fit.dx_fit/fit.f.a).*cos(fit.dx_fit/fit.f.b+fit.f.c)/cos(fit.f.c);
%cxt_pfit = polyval(pfit_p,dx_fit

% figure()
% scatter(dx_nond2_5loc_all,cxt_dxwidth_5loc_all,'k');
% hold on
% plot(dx_fit,cxt_fit,'LineWidth',2)
% hold on 
% errorbar(binmean2.dx_bincenter,binmean2.cxt_binmean,binmean2.error_binmean,'Color', 'm','LineWidth',2)
% hold off 
%%
save('/data1/bliu/data/cxt_dxwidth',"cxt_dxwidth","cxt_dxwidth_5loc",'binmean1','binmean2','fit')


