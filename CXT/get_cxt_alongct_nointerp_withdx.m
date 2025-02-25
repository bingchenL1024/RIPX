% Bingchen Liu Feb 13, 2024
% This code extract CXT along ct (exact, without consider magnitude of cxt)
% This code also extact cross-shore dx structure of cxt at each x=ct. Note
% have +/- of 1m or 2 meter


clear
close all

load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen_qced.mat')

load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
%%

for runnum = 1:120
    g = 9.81;
    dim = size(cell2mat(CXT_ALL(runnum)));
    SS = S(runnum);
    x = SS.X;
    x2=SS.X2;
    h = SS.h;
    Hs_raw = SS.Hs;
    Hs_interp= interp1(x2,Hs_raw,x);

    cxt = cell2mat(CXT_ALL(runnum));
    cxt_upbd=cell2mat(CXT_upbd_ALL(runnum));
    cxt_lowbd=cell2mat(CXT_lowbd_ALL(runnum));
    xb= SS.xb;
    x_nond = x./xb;
    G0_nond = abs(SS.curlF_std)./abs(max(SS.curlF_std));%G0 modification
    
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
        h_xloc = h(xind)+Hs_interp(xind);
        c = sqrt(g*h_xloc);
        ct = c.*t_itp;
        ct_all(:,i) = ct;
        cxt_atx = squeeze(cxt(xind,:,:));
        cxt_atx_uncert = squeeze(cxt_upbd(xind,:,:)-cxt_lowbd(xind,:,:));
        % if xind == 100
        %     keyboard
        % end 
        cxt_alongct(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct,'nearest');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% use BL close to ct creteria to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% select CXT 
        for j=1:length(t_itp)
            % if i ==50&& j == 3 &&runnum==120
            %     keyboard
            % end 
            [~,rowind_ct] = min(abs(dx(:,t_center+j)-ct(j))); % find index of row that is closest to x=ct
            %cxt_alongct(j,i) = cxt_atx(rowind_ct,t_center+j);
            cxt_alongct_pm1dx{j,i} = nan(1,3);
            cxt_alongct_pm2dx{j,i} = nan(1,5);
            if ~isnan(cxt_alongct(j,i))
                if j ==1 
                    cxt_alongct_pm1dx{j,i}(2:end) = cxt_atx(1:2,t_center+j);
                    cxt_alongct_pm2dx{j,i}(3:end) = cxt_atx(1:3,t_center+j);
                else
                    cxt_dxwidth_pm1_temp = cxt_atx(rowind_ct-1:min(rowind_ct+1,dim(2)),t_center+j);
                    cxt_dxwidth_pm2_temp = cxt_atx(rowind_ct-2:min(rowind_ct+2,dim(2)),t_center+j);
                    
                    cxt_alongct_pm1dx{j,i}(1:length(cxt_dxwidth_pm1_temp)) = cxt_dxwidth_pm1_temp;
                    cxt_alongct_pm2dx{j,i}(1:length(cxt_dxwidth_pm2_temp)) = cxt_dxwidth_pm2_temp;
                    
                    clear cxt_dxwidth_pm1_temp cxt_dxwidth_pm2_temp
                end 
            end
        end 
        
    end 
    cxt_alongct_ALL{runnum} = cxt_alongct;
    cxt_alongct_pm1dx_ALL{runnum} = cxt_alongct_pm1dx;   
    cxt_alongct_pm2dx_ALL{runnum} = cxt_alongct_pm2dx;   

    x_nond_All{runnum} =  x_nond(ind_good); 
    G0_nond_all{runnum} = G0_nond;

    
    % ind_nonzero  = find (cxt_alongct(1,:) ~=0);
    % cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
    %cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

%cxt_alongct_ALL = cxt_alongct_ALL';
%cxt_alongct_mean = cxt_alongct_mean';
%ind_good_All = ind_good_All';

dxwidth_cxtcoord_pm1= [-1,0,1];
dxwidth_cxtcoord_pm2= [-2,-1,0,1,2];

save('/data1/bliu/data/cxt_alongct_nointerp_dxwidth','cxt_alongct_ALL','cxt_alongct_pm1dx_ALL','cxt_alongct_pm2dx_ALL', ...
    'dxwidth_cxtcoord_pm1','dxwidth_cxtcoord_pm2','t_itp','x_nond_All')
