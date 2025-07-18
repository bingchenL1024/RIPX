% Bingchen Liu Sep 25, 2024
% This code plot cxt path with different phase velocity (c = 0.5~2 * sqrt(gh))

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
xb= SS.xb;
x_nond = x./xb;
G0 = SS.curlF_std;%G0 modification
cxt = cell2mat(CXT_ALL(runnum));
runpar = get_runpara(runnum);
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
ind_insz = find(x_nond>-1);
ind_bigG0 = find(abs(G0)>0.1*abs(max(G0)));%G0 modification
ind_good = intersect(ind_good,ind_insz);  % pick cxt inside the surf zone 
ind_good = intersect(ind_good,ind_bigG0); %G0 modification
ind_good_All{runnum} =ind_good; 

x_lag = 0:9; % in meters
t_lag = -50:50; %in seconds
    
[dt,dx] = meshgrid(t_lag,x_lag);

x_itp = 0:0.1:9; %interpolate into 0.1 m res
cxt_alongct= zeros(length(x_itp),length(ind_good));
t_all = zeros(length(x_itp),length(ind_good));

for i=1:5:length(ind_good) %cross shore dim
    xind = ind_good(i);
    h_xloc = h(xind);
    c = sqrt(g*h_xloc);
    c_gh_all(i) =c;
    t_gh_all(:,i) = x_itp./c_gh_all(i);
    cxt_atx = squeeze(cxt(xind,:,:));
    c_coef_all = 0.5+((1:31)-1).*0.1;
    for c_coef_ind = 1:31
        c_coef = 0.5+(c_coef_ind-1).*0.1;
        t_itp_temp= x_itp./(c.*c_coef);
        t_itp_allc(:,c_coef_ind) = t_itp_temp;
        cxt_alongct_temp = interp2(dt,dx,cxt_atx,t_itp_temp,x_itp);
        cxt_alongct_allc(:,c_coef_ind) = cxt_alongct_temp;
        cxt_var(c_coef_ind) = sum(cxt_alongct_temp.^2); %int (x^2)
        
        %debug =====================================
%         if i == 11
%         figure(1)
%         contourf(dt,dx,cxt_atx,20,'LineColor','none')
%         hold on 
%         plot(t_itp_temp,x_itp,'LineWidth',2,'color','k')
%         title(['var=',num2str(var(cxt_alongct_temp)),'   c=', num2str(c_coef)])
%         pause(1)
%         end 
        %==========================================
    end 

    [~,ind_maxvar] = max(cxt_var);
    ind_maxvar_all(i) = ind_maxvar;
    t_all(:,i) = t_itp_allc(:,ind_maxvar);
    %ct_all(:,i) = ct;
    cxt_alongct(:,i)= cxt_alongct_allc(:,ind_maxvar);
    c_coef_maxvar(i) = 0.5+(ind_maxvar-1).*0.1; 
    c_maxvar_all(i) = c.*c_coef_maxvar(i); 

    %=============================================================plot
    clf(figure(2))
    close all
    figure(2)
    sgtitle([runpar.wave,runpar.bath],'Fontsize',18)

    subplot(3,1,1)
    %contourf(dt,dx,cxt_atx,20,'LineColor','none')
    p1=pcolor(dt,dx, cxt_atx); %pcolor
    shading interp
    hold on 
    h1=plot(t_all(:,i),x_itp,'LineWidth',2,'color','r');
    hold on 
    h2=plot(t_gh_all(:,i),x_itp,'LineWidth',2,'color','b');
    ylim([0,9])
    xlim([-8,8])
    col=colorbar;
    cmocean('balance');
    %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
    caxis([-1,1])
    col.Label.String = 'Cross Correlation';
    legend([h1,h2],{'max var','sqrt(gh)'})
    xlabel('dt')
    ylabel('dx')
    title(['h=',num2str(h_xloc),'m','   coef=', num2str(c_coef_maxvar(i))...
        , ' G0/G0_{max}=',sprintf('%.2f\n',G0(xind)/max(G0))])
    niceplot_nobold_nomintick(16)
    grid on 
    
    subplot(3,1,2)
    scatter(c_coef_all,cxt_var,80, 'filled')
    xlabel('coef (c=coef*sqrt(gh))')
    ylabel('$\int (C_{XT})^2$','interpreter','latex');
    niceplot_nobold_nomintick(16)
    grid on 
    
    subplot(3,1,3)
    plot(x_itp,cxt_alongct_allc)
    hold on 
    p1 = plot(x_itp,cxt_alongct_allc(:,ind_maxvar),'o','LineWidth',2,'color','r');
    hold on 
    p2 = plot(x_itp,cxt_alongct_allc(:,6),'o','LineWidth',2,'Color','b');
    title(['h=',num2str(h_xloc),'m',' coef=',num2str(c_coef_maxvar(i)), ' G0/G0_{max}=',sprintf('%.2f\n',G0(xind)/max(G0))])
    legend([p1,p2],{'max var','sqrt(gh)'})
    xlabel('Space Lag (m)')
    ylabel('Cross Correlation')
    grid on 
    niceplot_nobold_nomintick(16)
     
    width= 6;
    height = 15;
    set(gcf,'Units','inches','Position',[0,0,width,height])
    set(gcf,'visible','off') 
    shg
    saveas(gcf,['/data1/bliu/figures/Cxt_coef/','run_',num2str(runnum),'xind_',num2str(xind),'.png'])
    %================================================= end of plot 
end 



end 