% Bingchen Liu Aug 23, 2024
% plot cxt(dt,dx)

%clear

%load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt.mat')
load ('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat')
load('/data1/bliu/data/SS_raw.mat')
load('/data1/bliu/data/cxt_ind_good')

%%
for runnum = 1:120
    everyxloc = 1;

    %%
    g = 9.81;
    dim = size(cell2mat(CXT_ALL(runnum)));
    SS = S(runnum);
    x = SS.X;
    h = SS.h;
    xb= SS.xb;
    x_nond = x./xb;
    cxt = cell2mat(CXT_ALL(runnum));

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
    % ind_good_All{runnum} =ind_good; 
    ind_good = ind_good_All{runnum};
    x_lag = 0:dim(2)-1; % in meters
    t_lag = -10:10; %in seconds

    [dt,dx] = meshgrid(t_lag,x_lag);
    %%
    for i=1:everyxloc:length(ind_good) %cross shore dim
        %xloc=  111;
        h_xloc = h(ind_good(i));
        c = sqrt(g*h_xloc);

        cxt_atx = squeeze(cxt(ind_good(i),:,:));

        ct = c.*t_lag;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%cxt(dx,dt)
        pcolorcen(dt,dx,cxt_atx);
        %pcolor(dt,dx, cxt_atx)
        hold on 
        plot(t_lag,ct,'LineWidth',2,'color','k')
        %shading interp


        col=colorbar;
        cmocean('balance');
        %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
        caxis([-1,1])
        col.Label.Interpreter = 'latex';
        col.Label.String = '$\tilde{C}_{XT}$';
        col.Label.FontSize = 22;
        col.Label.FontWeight = 'bold';
        xticks(-10:1:10)
        ylabel('$\Delta x$ (m)','Interpreter','latex')
        xlabel('$\Delta t$ (s)','Interpreter','latex')
        %title(['xlocation = ',num2str(x_nond(ind_good(i))),' dimensionless sz loc'])
        ylim([0,dim(2)])
        xlim([-10,10])
        niceplot(18);
    %     drawnow
    %     pause(1)
        width= 10;
        height = 10;
        set(gcf,'Units','inches','Position',[0,0,width,height]);
        set(gcf,'visible','off');
        disp([num2str(runnum),',',num2str(i)])
        saveas(gcf,['/data1/bliu/figures/RIPX_allrun/CXT_dtdx/nointerp_21dx10dt_1mres_BLnormalization/','run_',num2str(runnum),'_xloc_',num2str(ind_good(i)),'.png'])
        clf

    end 
end 
