% Bingchen Liu Nov 27, 2024
% This code save the data for snopshot of vorticity and curlF field for paper
% generate data used in 'plot_2d_vort_curlF'
runind = 14;
%runind = 44;
clearvars -except runind
t= 2000;


fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',runind);
load(fname,'curlF')
curlF_snap = squeeze(curlF(:,:,t))';
clear curlF


fname1=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_vort/RIPX_vort_%04d',runind);
load(fname1,'vort')
vort_snap = squeeze(vort(t,:,:))';
clear vort


head = 'data generated from "get_snap_vortcurlFfield.m"';
save('/data1/bliu/data/snap_vort_curlF_field_run14','vort_snap',"curlF_snap","head")
%% test 
% t= 2036;
% dim= size(curlF);
% x = 0:dim(1)-1;
% y = 0:dim(2)-1; %in seconds
% [x_grid,y_grid] = meshgrid(x,y);
% % for t = 2000:2100;
% figure(1)
% subplot(211)
% curlF_temp = curlF(:,:,t);
% pcolorcen(y_grid,x_grid,curlF_snap);
% col=colorbar;
% cmocean('balance');
% %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
% caxis([-0.5,0.5])
% %col.Limits = [-0.35,0.35];
% col.Label.Interpreter = 'latex';
% col.Label.String = '$\nabla \times F_{br} \, (s^{-2})$';
% col.Label.FontSize = 22;
% col.Label.FontWeight = 'bold';
% %col.TickLabels=[-0.35,-0.2,0,0.2,0.35];
% ylabel('x (m)','Interpreter','latex')
% xlabel('y (m)','Interpreter','latex')
% title(num2str(t))
% % width= 25;
% % height = 10;
% % set(gcf,'Units','inches','Position',[0,0,width,height]);
% % set(gcf,'visible','off');
% ylim([1,dim(1)])
% xlim([0,350])
% set(gca,'YDir','reverse')
% niceplot(22)
% 
% 
% subplot(212)
% vort_temp = vort(:,:,t);
% pcolorcen(y_grid,x_grid,vort_snap);
% col=colorbar;
% cmocean('balance');
% %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
% caxis([-0.1,0.1])
% col.Label.Interpreter = 'latex';
% col.Label.String = '$\nabla \times u \, (s^{-1})$';
% col.Label.FontSize = 22;
% col.Label.FontWeight = 'bold';
% ylabel('x (m)','Interpreter','latex')
% xlabel('y (m)','Interpreter','latex')
% %title([runpara.wave,runpara.bath])
% % width= 22;
% % height = 18;
% % set(gcf,'Units','inches','Position',[0,0,width,height]);
% % set(gcf,'visible','off');
% %ylim([dim(1)-100,dim(1)])
% ylim([1,dim(1)])
% xlim([0,350])
% set(gca,'YDir','reverse')
% niceplot(22)
% 
% % drawnow
% % pause(0.5)
% % end 