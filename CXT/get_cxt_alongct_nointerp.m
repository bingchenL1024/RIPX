% Bingchen Liu Nov 5, 2024
% This code extract CXT along ct 
% IMPORTANT: it does NOT include interp of cxt in noninterger t but only
% extract with discrete t, i.e. t=0, 1, 2, 3...


clear
close all

%load ('/data1/bliu/data/raw/CXT_ALL_Mark.mat')
%load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')

load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
%%

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



    
    % cxt_min = zeros(dim(1),1);
    % cxt_max = zeros(dim(1),1);
    % for i = 1:dim(1)
    %     slice = cxt(i,:,:);
    %     cxt_max(i) = max(slice(:));
    %     cxt_min(i) = min(slice(:));
    % end 
    % ind_good_temp = find(cxt_max<=1&cxt_min>=-1&cxt_max>0&cxt_min<0);
    % ind_good_diff = diff(ind_good_temp);
    % ind_good = ind_good_temp(find(ind_good_diff==1));
    % ind_insz = find(x_nond>-1);
    % ind_good = intersect(ind_good,ind_insz);  % pick cxt inside the surf zone 
    % ind_good_All{runnum} = ind_good; 
    
    ind_good = ind_good_All{runnum};
    
    x_lag = 0:20; % in meters
    t_lag = -10:10; %in seconds
    t_center = (length(t_lag)-1)/2;
    [dt,dx] = meshgrid(t_lag,x_lag);
    t_itp = 0:1:10;
    cxt_alongct= zeros(length(t_itp),length(ind_good));
    ct_all = zeros(length(t_itp),length(ind_good));
    
    for i=1:length(ind_good) %cross shore dim
        xind = ind_good(i);
        h_xloc = h(xind);
        c = sqrt(g*h_xloc);
        ct = c.*t_itp;
        ct_all(:,i) = ct;
        cxt_atx = squeeze(cxt(xind,:,:));
        cxt_atx_uncert = squeeze(cxt_upbd(xind,:,:)-cxt_lowbd(xind,:,:));
        if xind == 100
            keyboard
        end 
        cxt_interp(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct);
        threshold_cxt = ceil(c*0.5); % threshold on how close the data needs to be close to x=ct (include c variation since 1s cover diff distance),otherwise, use interpolation.
        range_maxcxt = ceil(c*1); %the range where we look for the max CXT around x= ct
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% use BL close to ct creteria to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% select CXT 
        %[val_max,ind_max] = max(cxt_atx);
        for j=1:length(t_itp)
            % if i ==50&& j == 3 &&runnum==120
            %     keyboard
            % end 
            [~,rowind_ct] = min(abs(dx(:,t_center+j)-ct(j))); % find index of row that is closest to x=ct
            cxt_aroundct = cxt_atx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get cxt at given time and cross-shore loc around x=ct
            dx_aroundct =dx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get the same info as above but for dx
            [val_max,ind_max] = max(cxt_aroundct);
            if abs(dx_aroundct(ind_max)-ct(j))<threshold_cxt %if the max value is very close (<1) to x=ct 
                cxt_alongct(j,i) = val_max;
                x_cxt_alongct(j,i) = dx_aroundct(ind_max);
                if_interp(j,i)=0;
            else 
                cxt_alongct(j,i)=cxt_interp(j,i);
                if ~isnan(cxt_interp(j,i))
                    x_cxt_alongct(j,i) = ct(j);
                else 
                    x_cxt_alongct(j,i) = NaN;
                end 
                if_interp(j,i) = 1;
            end 
            clear cxt_aroundct dx_aroundct val_max ind_max
        end 
        
    end 
    if_interp_All{runnum} = if_interp;
    clear if_interp
    cxt_alongct_ALL{runnum} = cxt_alongct;
    x_cxt_alongct_All{runnum} = x_cxt_alongct;
    x_nond_All{runnum} =  x_nond(ind_good); 
    G0_nond_all{runnum} = G0_nond;

    
    % ind_nonzero  = find (cxt_alongct(1,:) ~=0);
    % cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
    %cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

%cxt_alongct_ALL = cxt_alongct_ALL';
%cxt_alongct_mean = cxt_alongct_mean';
%ind_good_All = ind_good_All';

save('/data1/bliu/data/cxt_alongct_nointerp','cxt_alongct_ALL','if_interp_All','x_cxt_alongct_All','t_itp','x_nond_All')
