%Bingchen Liu July 26, 2024
% This code extract cxt along x=ct line and it loops through different run 
% Note it uses space x as x axis 
% Note it interpolate x into 0.1 m resolution
% Note it only extract cxt inside of the surf zone 

%Major modification Sep 18
%*search for group velocity c st the variance of cxt along x=ct is max

clear
close all

%load ('/data1/bliu/data/raw/CXT_ALL_Mark.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_Falk.mat')
%load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')

load('/data1/bliu/data/SS_raw.mat')

G0_lim = 0.15;

for runnum = 1:120
g = 9.81;
dim = size(cell2mat(CXT_ALL(runnum)));
SS = S(runnum);
x = SS.X;
h = SS.h;
xb= SS.xb;
x_nond = x./xb;
G0 = SS.curlF_std;%G0 modification
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
ind_bigG0 = find(abs(G0)>G0_lim*abs(max(G0)));%G0 modification
ind_insz = find(x_nond>-1);
ind_good = intersect(ind_good,ind_insz);  % pick cxt inside the surf zone 
ind_good = intersect(ind_good,ind_bigG0); %G0 modification
ind_good_All{runnum} =ind_good; 

x_lag = 0:9; % in meters
t_lag = -50:50; %in seconds
    
[dt,dx] = meshgrid(t_lag,x_lag);

x_itp = 0:0.1:9; %interpolate into 0.1 m res
cxt_alongct= zeros(length(x_itp),length(ind_good));
t_all = zeros(length(x_itp),length(ind_good));

for i=1:length(ind_good) %cross shore dim
    xind = ind_good(i);
    h_xloc = h(xind);
    c = sqrt(g*h_xloc);
    c_gh_all(i) =c;
    t_all(:,i) = x_itp./c_gh_all(i);
    cxt_atx = squeeze(cxt(xind,:,:));
    cxt_alongct_atx = interp2(dt,dx,cxt_atx,t_all(:,i)',x_itp);
    cxt_alongct(:,i) = cxt_alongct_atx;
    %ct_all(:,i) = ct;
end 

cxt_alongct_ALL{runnum} = cxt_alongct;
x_nond_All{runnum} =  x_nond(ind_good); 
t_itp_all{runnum} = t_all;
% ind_nonzero  = find (cxt_alongct(1,:) ~=0);
% cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
% cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

cxt_alongct_x_ALL = cxt_alongct_ALL';
%cxt_alongct_mean = cxt_alongct_mean';
x_nond_All = x_nond_All';
ind_good_All = ind_good_All';

save('/data1/bliu/data/cxt_alongct_x','cxt_alongct_x_ALL','t_itp_all','x_itp','x_nond_All','ind_good_All')

% for i = 1:120
% plot(t_itp,cxt_alongct_mean(:,i))
% hold on 
% end 
% hold off