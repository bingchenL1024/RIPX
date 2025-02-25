%Bingchen Liu July 26, 2024
% This code extract cxt along x=ct line and it loops through different run 
% Note it uses time as x axis


clear
close all

%load ('/data1/bliu/data/raw/CXT_ALL_Mark.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_Falk.mat')
load('/data1/bliu/data/SS_raw.mat')

for runnum = 1:120
g = 9.81;
dim = size(cell2mat(CXT_ALL(runnum)));
SS = S(runnum);
x = SS.X;
h = SS.h;
cxt = cell2mat(CXT_ALL(runnum));

cxt_min = zeros(dim(1),1);
cxt_max = zeros(dim(1),1);
for i = 1:dim(1)
    slice = cxt(i,:,:);
    cxt_max(i) = max(slice(:));
    cxt_min(i) = min(slice(:));
end 
ind_good_temp = find(cxt_max<=1&cxt_min>=-1&cxt_max>0&cxt_min<0);
ind_good_diff = diff(ind_good_temp);
ind_good = ind_good_temp(find(ind_good_diff==1));
ind_good_All{runnum} = ind_good; 


x_lag = 0:9; % in meters
t_lag = -50:50; %in seconds
    
[dt,dx] = meshgrid(t_lag,x_lag);

t_itp = 0:0.1:10;
cxt_alongct= zeros(length(t_itp),length(ind_good));
ct_all = zeros(length(t_itp),length(ind_good));

for i=1:length(ind_good) %cross shore dim
    xind = ind_good(i);
    h_xloc = h(xind);
    c = sqrt(g*h_xloc);
    ct = c.*t_itp;
    ct_all(:,i) = ct;
    cxt_atx = squeeze(cxt(xind,:,:));
    cxt_alongct(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct);
    
end 

cxt_alongct_ALL{runnum} = cxt_alongct;

ind_nonzero  = find (cxt_alongct(1,:) ~=0);
cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
%cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

cxt_alongct_t_ALL = cxt_alongct_ALL';
cxt_alongct_mean = cxt_alongct_mean';
ind_good_All = ind_good_All';

save('/data1/bliu/data/cxt_alongct_t','cxt_alongct_t_ALL','cxt_alongct_mean','t_itp','ind_good_All')

% for i = 1:120
% plot(t_itp,cxt_alongct_mean(:,i))
% hold on 
% end 
% hold off