% Bingchen Liu Feb 19, 2025
% This code extract CXT along ct 
% IMPORTANT: it does NOT include interp of cxt in noninterger t but only
% extract with discrete t, i.e. t=0, 1, 2, 3...
% Note it's different from the original nointerp max code that we don't do
% something different when dx is not close enough to c*dt but deal with
% that fomr RMSE creteria 
% include RMSE of dx-c*dt
% include dx width 


clear
close all

load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen_qced.mat')
load('/data1/bliu/data/runnum_72run')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
%%

CXT_threshold = 0.025;

for runnum = 1:120
    g = 9.81;
    dim = size(cell2mat(CXT_ALL(runnum)));
    SS = S(runnum);
    x = SS.X;
    h = SS.h;
    cxt = cell2mat(CXT_ALL(runnum));
    cxt_upbd=cell2mat(CXT_upbd_ALL(runnum));
    cxt_lowbd=cell2mat(CXT_lowbd_ALL(runnum));
    xb= SS.xb;
    x_nond = x./xb;
    G0_nond = abs(SS.curlF_std)./abs(max(SS.curlF_std));%G0 modification
    
    ind_good = ind_good_All{runnum};
    
    x_lag = 0:20; % in meters HARD CODED
    t_lag = -10:10; %in seconds  HARD CODED
    t_center = (length(t_lag)-1)/2;
    [dt,dx] = meshgrid(t_lag,x_lag);
    t_itp = 0:1:5; % HARD CODED
    cxt_alongct= NaN(length(t_itp),length(ind_good));
    ct_all = NaN(length(t_itp),length(ind_good));
    dist2ct = NaN(length(t_itp),length(ind_good));
    dx_max = NaN(length(t_itp),length(ind_good));
    dx_ct = NaN(length(t_itp),length(ind_good));
    for i=1:length(ind_good) %cross shore dim
        xind = ind_good(i);
        h_xloc = h(xind);
        c = sqrt(g*h_xloc);
        ct = c.*t_itp;
        ct_all(:,i) = ct;
        cxt_atx = squeeze(cxt(xind,:,:));
        cxt_atx_uncert = squeeze(cxt_upbd(xind,:,:)-cxt_lowbd(xind,:,:));




        cxt_interp(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct);
        threshold_cxt = ceil(c*0.5); % threshold on how close the data needs to be close to x=ct (include c variation since 1s cover diff distance),otherwise, use interpolation.
        range_maxcxt = ceil(c*5); %the range where we look for the max CXT around x= ct
    

        for j=1:length(t_itp)
        % begug ================================
        % if runnum==98&&i ==35 && j ==6
        %     plot_CXT_dxdt_function(cxt_atx,c)
        %     keyboard
        % 
        % end 
        % begug ================================
         
            [~,rowind_ct] = min(abs(dx(:,t_center+j)-ct(j))); % find index of row that is closest to x=ct
            cxt_aroundct = cxt_atx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get cxt at given time and cross-shore loc around x=ct
            dx_aroundct =dx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get the same info as above but for dx
            [val_max,ind_max] = max(cxt_aroundct);
            [~,ind_max_fulldxdomain] = min(abs(cxt_atx(:,t_center+j)-val_max)); % find the index of the CXT we pick relative to whole cxt_atx domain (given t)
            
            cxt_alongct_pm1dx{j,i} = nan(1,3);
            cxt_alongct_pm2dx{j,i} = nan(1,5);
            if ct(j)<x_lag(end) % make sure the pt is in the domain
                if j ==1 && val_max>CXT_threshold% at 0 tlag 
                    cxt_alongct(j,i) = val_max;
                    if ~isnan(cxt_alongct(j,i)) % only record location when cxt is not NaN
                        dist2ct(j,i) = dx_aroundct(ind_max)-ct(j);
                        dx_max(j,i) = dx_aroundct(ind_max);
                        dx_ct(j,i) = ct(j);
                    end 
                    cxt_alongct_pm1dx{j,i}(2:end) = cxt_atx(1:2,t_center+j); % 2 is center in dx coord
                    cxt_alongct_pm2dx{j,i}(3:end) = cxt_atx(1:3,t_center+j);% 3 is center in dx coord
                elseif val_max>CXT_threshold
                    cxt_alongct(j,i) = val_max;
                    if ~isnan(cxt_alongct(j,i)) % only record location when cxt is not NaN
                        dist2ct(j,i) = dx_aroundct(ind_max)-ct(j);
                        dx_max(j,i) = dx_aroundct(ind_max);
                        dx_ct(j,i) = ct(j);
                    end %test 
                    cxt_dxwidth_pm1_temp = nan(3,1);
                    cxt_dxwidth_pm2_temp = nan(5,1);

                    cxt_dxwidth_pm1_temp = cxt_atx(ind_max_fulldxdomain-1:min(ind_max_fulldxdomain+1,dim(2)),t_center+j);
                    cxt_dxwidth_pm2_temp = cxt_atx(ind_max_fulldxdomain-2:min(ind_max_fulldxdomain+2,dim(2)),t_center+j);
                    cxt_alongct_pm1dx{j,i}(1:length(cxt_dxwidth_pm1_temp)) = cxt_dxwidth_pm1_temp;
                    cxt_alongct_pm2dx{j,i}(1:length(cxt_dxwidth_pm2_temp)) = cxt_dxwidth_pm2_temp;
                    
                    clear cxt_dxwidth_pm1_temp cxt_dxwidth_pm2_temp
                end 
            end 
            clear cxt_aroundct dx_aroundct val_max ind_max
        end 
    end 
    cxt_alongct_ALL{runnum} = cxt_alongct;
    
    cxt_alongct_pm1dx_ALL{runnum} = cxt_alongct_pm1dx;   
    cxt_alongct_pm2dx_ALL{runnum} = cxt_alongct_pm2dx;
    
    dist2ct_All{runnum} = dist2ct;
    %MSE_dist(runnum) = var(dist2ct(:),'omitmissing');
    dx_max_All{runnum} = dx_max; %location of max CXT
    dx_ct_All{runnum} = dx_ct; % location of ct 
    ind_nonan = ~isnan(dx_max)&~isnan(dx_ct);
    MSE_dist(runnum) = immse(dx_max(ind_nonan),dx_ct(ind_nonan));
    clear ind_nonan

    x_nond_All{runnum} =  x_nond(ind_good); 
    G0_nond_all{runnum} = G0_nond;

    
    % ind_nonzero  = find (cxt_alongct(1,:) ~=0);
    % cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
    %cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

dxwidth_cxtcoord_pm1= [-1,0,1];
dxwidth_cxtcoord_pm2= [-2,-1,0,1,2];
RMSE_dist = sqrt(mean(MSE_dist(goodrunnum)));



%cxt_alongct_ALL = cxt_alongct_ALL';
%cxt_alongct_mean = cxt_alongct_mean';
%ind_good_All = ind_good_All';

save('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth','cxt_alongct_ALL','cxt_alongct_pm1dx_ALL','cxt_alongct_pm2dx_ALL', ...
    'dist2ct_All','t_itp','x_nond_All','RMSE_dist','dxwidth_cxtcoord_pm1','dxwidth_cxtcoord_pm2')
