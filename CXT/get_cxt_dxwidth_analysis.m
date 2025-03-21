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

dx_width = -2:1:2;

for ind = 1:24
    runnum = indbath.slp2(ind);
    a_decotscale = fitpara.slp2(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp2{ind,1}.dx_nond1{1,xloc} = dx_width/(c_modelh{runnum}(xloc)*fitpara.slp2(ind).a(xloc));
        cxt_dxwidth.slp2{ind,1}.dx_nond2{1,xloc} = dx_width/(c_modelh{runnum}(xloc));
        cxt_dxwidth.slp2{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp2{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

    end 
    cxt_dxwidth.slp2{ind,1}.x_nond = x_nond_All{runnum};
end 


for ind = 1:24
    runnum = indbath.slp3(ind);
    a_decotscale = fitpara.slp3(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp3{ind,1}.dx_nond1{1,xloc} = dx_width/(c_modelh{runnum}(xloc)*fitpara.slp3(ind).a(xloc));
        cxt_dxwidth.slp3{ind,1}.dx_nond2{1,xloc} = dx_width/(c_modelh{runnum}(xloc));
        cxt_dxwidth.slp3{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp3{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

    end 
    cxt_dxwidth.slp3{ind,1}.x_nond = x_nond_All{runnum};

end 

for ind = 1:24
    runnum = indbath.slp4(ind);
    a_decotscale = fitpara.slp4(ind).a;
    for xloc = 1:length(a_decotscale)
        cxt_dxwidth.slp4{ind,1}.dx_nond1{1,xloc} = dx_width/(c_modelh{runnum}(xloc)*fitpara.slp4(ind).a(xloc));
        cxt_dxwidth.slp4{ind,1}.dx_nond2{1,xloc} = dx_width/(c_modelh{runnum}(xloc));
        cxt_dxwidth.slp4{ind,1}.cxt_nond{1,xloc} = cxt_pm2dx_nond_ALL{runnum}{1,xloc};
        cxt_dxwidth.slp4{ind,1}.cxt{1,xloc} = cxt_pm2dx_ALL{runnum}{1,xloc};

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
        cxt_dxwidth_5loc.slp2(i).cxt(j,:) = cxt_dxwidth.slp2(i).cxt{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp2(i).cxt_nond(j,:) = cxt_dxwidth.slp2(i).cxt_nond{ind_5loc_temp(j)};

    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(cxt_dxwidth.slp3(i).x_nond - loc5(j)));
    
        cxt_dxwidth_5loc.slp3(i).x_nond(j) = cxt_dxwidth.slp3(i).x_nond(ind_5loc_temp(j));
        cxt_dxwidth_5loc.slp3(i).dx_nond1(j,:) = cxt_dxwidth.slp3(i).dx_nond1{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp3(i).dx_nond2(j,:) = cxt_dxwidth.slp3(i).dx_nond2{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp3(i).cxt(j,:) = cxt_dxwidth.slp3(i).cxt{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp3(i).cxt_nond(j,:) = cxt_dxwidth.slp3(i).cxt_nond{ind_5loc_temp(j)};

    end 

    for j = 1:5
        [~,ind_5loc_temp(j)] = min(abs(cxt_dxwidth.slp4(i).x_nond - loc5(j)));
    
        cxt_dxwidth_5loc.slp4(i).x_nond(j) = cxt_dxwidth.slp4(i).x_nond(ind_5loc_temp(j));
        cxt_dxwidth_5loc.slp4(i).dx_nond1(j,:) = cxt_dxwidth.slp4(i).dx_nond1{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp4(i).dx_nond2(j,:) = cxt_dxwidth.slp4(i).dx_nond2{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp4(i).cxt(j,:) = cxt_dxwidth.slp4(i).cxt{ind_5loc_temp(j)};
        cxt_dxwidth_5loc.slp4(i).cxt_nond(j,:) = cxt_dxwidth.slp4(i).cxt_nond{ind_5loc_temp(j)};

    end 



end 

%% take avg of dx structure at 5 loc 


dxlag = 3:5;
cxt_dxwidth_all = [];
for ind = 1:24
        cxt_temp = [cxt_dxwidth_5loc.slp2(ind).cxt(:,dxlag);cxt_dxwidth_5loc.slp3(ind).cxt(:,dxlag);cxt_dxwidth_5loc.slp4(ind).cxt(:,dxlag)];
        cxt_dxwidth_all = [cxt_dxwidth_all;cxt_temp];
     
end 

save('/data1/bliu/data/cxt_dxwidth',"cxt_dxwidth","cxt_dxwidth_5loc",'cxt_dxwidth_all')


