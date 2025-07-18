% Bingchen Liu Feb 19, 2025
% This code extract CXT along ct 
% IMPORTANT: it does NOT include interp of cxt in noninterger t but only
% extract with discrete t, i.e. t=0, 1, 2, 3...
% Note it's different from the original nointerp max code that we don't do
% something different when dx is not close enough to c*dt but deal with
% that fomr RMSE creteria 
% include RMSE of dx-c*dt
% include dx width 

% April 28, 2025 update
% add interp using fit velocity c and interp around c*dt 

clear
close all

load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen_qced.mat')
load('/data1/bliu/data/runnum_72run')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')
load('/data1/bliu/data/ind_of_diff_bath.mat')

%%

CXT_threshold = 0; %0.025 previously 
t_itp = 0:1:5; % HARD CODED

Stats.N_dof_RMSE_dist = 0;%number of DOF for the dist RMSE calculation

for runnum = 1:120
    g = 9.81;
    dim = size(cell2mat(CXT_ALL(runnum)));
    SS = S(runnum);
    x = SS.X;
    x2= SS.X2;
    h = SS.h;
    hb = SS.hb;
    Hs = SS.Hs;
    ds_b= SS.sigma_b;
    Hs_interp = interp1(x2,Hs,x);
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
    cxt_alongct= NaN(length(t_itp),length(ind_good));
    cxt_alongct_interp= NaN(length(t_itp),length(ind_good));
    ct_all = NaN(length(t_itp),length(ind_good));
    dist2ct = NaN(length(t_itp),length(ind_good));
    dx_max = NaN(length(t_itp),length(ind_good));
    dx_ct = NaN(length(t_itp),length(ind_good));
    for i=1:length(ind_good) %cross shore dim
        xind = ind_good(i);
        h_xloc = h(xind);
        dirspr_b{runnum}(1,i) = ds_b;
        Hs_xloc= Hs_interp(xind);
        c = sqrt(g*h_xloc);
        c_nonlin_temp = sqrt(g*(h_xloc+Hs_xloc));
        c_modelh{runnum}(1,i) = c;
        c_modelhHs{runnum}(1,i) = c_nonlin_temp; 
        Hs_h{runnum}(1,i) = Hs_xloc/h_xloc; 
        h_hb_nond{runnum}(1,i) = h_xloc/hb;
        Hs_h_1psqrt{runnum}(1,i) = sqrt(1+Hs_xloc/h_xloc);
        

        ct = c.*t_itp;
        ct_all(:,i) = ct;
        cxt_atx = squeeze(cxt(xind,:,:));
        cxt_atx_uncert = squeeze(cxt_upbd(xind,:,:)-cxt_lowbd(xind,:,:));




        cxt_interp(:,i)= interp2(dt,dx,cxt_atx,t_itp,ct);
        threshold_cxt = ceil(c*0.5); % threshold on how close the data needs to be close to x=ct (include c variation since 1s cover diff distance),otherwise, use interpolation.
        range_maxcxt = ceil(c*5); %the range where we look for the max CXT around x= ct -- t = 5 so basically everywhere
    

        for j=1:length(t_itp)

         
            [~,rowind_ct] = min(abs(dx(:,t_center+j)-ct(j))); % find index of row that is closest to x=ct
            cxt_aroundct = cxt_atx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get cxt at given time and cross-shore loc around x=ct
            dx_aroundct =dx(max(rowind_ct-range_maxcxt,1):min(rowind_ct+range_maxcxt,max(x_lag)+1),t_center+j); %get the same info as above but for dx
            [val_max,ind_max] = max(cxt_aroundct);
            [~,ind_max_fulldxdomain] = min(abs(cxt_atx(:,t_center+j)-val_max)); % find the index of the CXT we pick relative to whole cxt_atx domain (given t)
            
            cxt_alongct_pm1dx{j,i} = nan(1,3);
            %cxt_alongct_pm2dx{j,i} = nan(5,1);
            cxt_alongct_pm2dx{j,i} = nan(9,1); % BL change dx to 0:4

            cxt_alongct_pm1dx_nond{j,i}= nan(1,3);
            % cxt_alongct_pm2dx_nond{j,i}= nan(1,5);%original
            cxt_alongct_pm2dx_nond{j,i}= nan(9,1); % BL change dx to 0:4


            if ct(j)<x_lag(end) % make sure the pt is in the domain
                if val_max>CXT_threshold %only record if CXT exceed threshold
                    cxt_alongct(j,i) = val_max;
                    if ~isnan(cxt_alongct(j,i)) % only record location when cxt is not NaN
                        dist2ct(j,i) = dx_aroundct(ind_max)-ct(j);
                        dx_max(j,i) = dx_aroundct(ind_max);
                        dx_ct(j,i) = ct(j);
                    end  
                    cxt_dxwidth_pm1_temp = nan(3,1);
                    % cxt_dxwidth_pm2_temp = nan(5,1); %original 
                    cxt_dxwidth_pm2_temp = nan(9,1); % BL change dx to 0:4
                    
                    
                    cxt_dxwidth_pm1_temp = extract_pm1(cxt_atx(:,t_center+j),ind_max_fulldxdomain);
                    cxt_dxwidth_pm2_temp = extract_pm2(cxt_atx(:,t_center+j),ind_max_fulldxdomain);
                    cxt_alongct_pm1dx{j,i} = cxt_dxwidth_pm1_temp;
                    cxt_alongct_pm2dx{j,i} = cxt_dxwidth_pm2_temp;


                    cxt_alongct_pm1dx_nond{j,i} = cxt_dxwidth_pm1_temp./cxt_dxwidth_pm1_temp(2);
                    % cxt_alongct_pm2dx_nond{j,i} = cxt_dxwidth_pm2_temp./cxt_dxwidth_pm2_temp(3);  %original                                  
                    cxt_alongct_pm2dx_nond{j,i} = cxt_dxwidth_pm2_temp./cxt_dxwidth_pm2_temp(5);   % BL change dx to 0:4                                 

                    clear cxt_dxwidth_pm1_temp cxt_dxwidth_pm2_temp
                end 
            end
            clear cxt_aroundct dx_aroundct val_max ind_max
        end
        
        cfit_ind = ~isnan(dx_max(:,i));
        dx_cfit = dx_max(cfit_ind,i);
        dt_cfit = t_itp(cfit_ind)';
        if length(dx_cfit)>=3 %only fit if we have more than 3 pt 
            c_fit_temp = dt_cfit\dx_cfit;
        else
            c_fit_temp= NaN;
        end 
        % begug ================================
        % if runnum==2&&i ==2 %&& j ==6
        %     plot_CXT_dxdt_function(cxt_atx,c)
        %     dt_cfit
        %     dx_cfit
        %     keyboard
        % 
        % end 
        % begug ================================
        c_fit{runnum}(1,i) = c_fit_temp;
        ct_itp = c_fit_temp.*t_itp;

        % April 28, 2025 update, use fit c to obtain cxt_alongct_interp
        % with interp 
       
        for j=1:length(t_itp)
            t_ind = t_center+j;
            cxt_alongct_interp(j,i)= interp1(x_lag,cxt_atx(:,t_ind),ct_itp(j));
        end 
            % if i==2
            %     keyboard
            % end 
    end 
    cxt_alongct_ALL{runnum} = cxt_alongct;
    cxt_alongct_itp_ALL{runnum} = cxt_alongct_interp;

    
    cxt_pm1dx_ALL{runnum} = cxt_alongct_pm1dx;   
    cxt_pm2dx_ALL{runnum} = cxt_alongct_pm2dx;

    cxt_pm1dx_nond_ALL{runnum} = cxt_alongct_pm1dx_nond;
    cxt_pm2dx_nond_ALL{runnum} = cxt_alongct_pm2dx_nond;


    
    dist2ct_All{runnum} = dist2ct;
    %MSE_dist(runnum) = var(dist2ct(:),'omitmissing');
    dx_max_All{runnum} = dx_max; %location of max CXT
    dx_ct_All{runnum} = dx_ct; % location of ct 
    ind_nonan = ~isnan(dx_max)&~isnan(dx_ct);
    Stats.N_dof_RMSE_dist = Stats.N_dof_RMSE_dist + sum(ind_nonan(:));
    MSE_dist(runnum) = immse(dx_max(ind_nonan),dx_ct(ind_nonan));
    clear ind_nonan

    x_nond_All{runnum} =  x_nond(ind_good); 
    G0_nond_all{runnum} = G0_nond;

    
    % ind_nonzero  = find (cxt_alongct(1,:) ~=0);
    % cxt_alongct_mean{runnum} = mean(cxt_alongct(:,ind_nonzero),2,'omitnan');
    %cxt_alongct_25_75 = prctile(cxt_alongct,[25 75],2);

end 

dxwidth_cxtcoord_pm1= [-1,0,1];
dxwidth_cxtcoord_pm2= -4:1:4;
Stats.RMSE_dist = sqrt(mean(MSE_dist(goodrunnum)));

%% extract all parameters at 5 surfzone location
loc5 = linspace(-0.75,-0.25,5);

for i = 1:24
    for j = 1:5
        runnum = indbath.slp2(i);
        [~,ind_5loc_temp(j)] = min(abs(x_nond_All{runnum} - loc5(j)));
        c_phase_5loc.slp2(i).dirspr_b(j) = dirspr_b{runnum}(ind_5loc_temp(j));

        c_phase_5loc.slp2(i).c_fit(j) = c_fit{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).c_modelh(j) = c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).c_modelhHs(j) = c_modelhHs{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).c_nond(j) = c_fit{runnum}(ind_5loc_temp(j))/c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).Hs_h(j) = Hs_h{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).Hs_h_sqrt(j) = Hs_h_1psqrt{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp2(i).h_hb_nond(j) = h_hb_nond{runnum}(ind_5loc_temp(j));

    end 

    for j = 1:5
        runnum = indbath.slp3(i);
        [~,ind_5loc_temp(j)] = min(abs(x_nond_All{runnum} - loc5(j)));
        c_phase_5loc.slp3(i).dirspr_b(j) = dirspr_b{runnum}(ind_5loc_temp(j));

        c_phase_5loc.slp3(i).c_fit(j) = c_fit{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).c_modelh(j) = c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).c_modelhHs(j) = c_modelhHs{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).c_nond(j) = c_fit{runnum}(ind_5loc_temp(j))/c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).Hs_h(j) = Hs_h{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).Hs_h_sqrt(j) = Hs_h_1psqrt{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp3(i).h_hb_nond(j) = h_hb_nond{runnum}(ind_5loc_temp(j));

    end 

    for j = 1:5
        runnum = indbath.slp4(i);
        [~,ind_5loc_temp(j)] = min(abs(x_nond_All{runnum} - loc5(j)));
        c_phase_5loc.slp4(i).dirspr_b(j) = dirspr_b{runnum}(ind_5loc_temp(j));

        c_phase_5loc.slp4(i).c_fit(j) = c_fit{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).c_modelh(j) = c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).c_modelhHs(j) = c_modelhHs{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).c_nond(j) = c_fit{runnum}(ind_5loc_temp(j))/c_modelh{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).Hs_h(j) = Hs_h{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).Hs_h_sqrt(j) = Hs_h_1psqrt{runnum}(ind_5loc_temp(j));
        c_phase_5loc.slp4(i).h_hb_nond(j) = h_hb_nond{runnum}(ind_5loc_temp(j));

    end 

end 


c_fit_tot=[];
c_modelh_tot=[];
for i= 1:24
    c_fit_tot= [c_fit_tot, c_phase_5loc.slp2(i).c_fit, c_phase_5loc.slp3(i).c_fit, c_phase_5loc.slp4(i).c_fit];
    c_modelh_tot= [c_modelh_tot, c_phase_5loc.slp2(i).c_modelh, c_phase_5loc.slp3(i).c_modelh, c_phase_5loc.slp4(i).c_modelh];

end 

nonan = ~isnan(c_fit_tot)&~isnan(c_modelh_tot);
Stats.RMSE_vel = rmse(c_fit_tot(nonan),c_modelh_tot(nonan));
Stats.c_corr= corrcoef(c_fit_tot(nonan),c_modelh_tot(nonan)).^2;
%cxt_alongct_ALL = cxt_alongct_ALL';
%cxt_alongct_mean = cxt_alongct_mean';
%ind_good_All = ind_good_All';

save('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth','cxt_alongct_ALL','cxt_pm1dx_ALL','cxt_pm2dx_ALL', ...
    'dist2ct_All','t_itp','x_nond_All','Stats','dxwidth_cxtcoord_pm1','dxwidth_cxtcoord_pm2', ...
    'cxt_pm1dx_nond_ALL','cxt_pm2dx_nond_ALL','c_modelh','c_fit','c_modelhHs','c_phase_5loc', ...
    'dx_max_All')
