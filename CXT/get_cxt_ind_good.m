% Bingchen Liu Nov 12, 2024
% This code get the cross-shore index x that is good based on creteria of
% cxt
% ******cxt(x) must be inside of the surf zone 
% ******cxt(x) needs to be between 0 and 1 and include negative and positive
% value
% ******correspondsing G0(x) of cxt(x) need to be biegger than G0_lim of G0max

clear
%load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/runnum_72run')

%%
G0_lim = 0.15;
test_num = 0;

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
    [ind_5loc] = get_10locs(SS);

    cxt_min = zeros(dim(1),1);
    cxt_max = zeros(dim(1),1);
    for i = 1:dim(1) %cross-shore x locations, every 1m 
        slice = cxt(i,:,:);
        cxt_max(i) = max(slice(:));
        cxt_min(i) = min(slice(:));
    end 
    % if runnum == 3
    %     keyboard 
    % end 
    ind_good_temp = find(cxt_max<=1&cxt_min>=-1&cxt_max>0&cxt_min<0);
    ind_bigG0 = find(abs(G0)>G0_lim*abs(max(G0)));%G0 modification
    ind_good_diff = diff(ind_good_temp);
    ind_good = ind_good_temp(find(ind_good_diff==1));
    ind_insz = find(x_nond>-1);
    ind_good = intersect(ind_good,ind_insz);  % pick cxt inside the surf zone 
    ind_good = intersect(ind_good,ind_bigG0); %G0 modification


    if ~all(ismember(ind_5loc,ind_good)) & ismember(runnum,goodrunnum)
        test_num = test_num+1;
        badrunnum(test_num) = runnum;
        test_diff(test_num)= ind_5loc(end)-ind_good(end);
    end 
    
    ind_good_All{runnum} =ind_good; 

end 



save('/data1/bliu/data/cxt_ind_good',"ind_good_All")