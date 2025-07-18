% Bingchen Liu, Nov 7, 2024
% this code produce the curlFbr video like the one emma produced 

clear all 
%for runind = 3:120
    runind=2;
    tic
    clearvars -except runind
    fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',runind);
    load(fname,'curlF')
    
    fname1=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_vort/RIPX_vort_%04d',runind);
    load(fname1,'vort')
    vort = permute(vort,[2 3 1]);
    %%
    close all
    runpara= get_runpara(runind);
    dim= size(curlF);
    x = 0:dim(1)-1;
    y = 0:dim(2)-1; %in seconds
    [x_grid,y_grid] = meshgrid(x,y);
    
    outputVideo = VideoWriter(['/data1/bliu/figures/Video/largedomain/run_',num2str(runind),runpara.wave,runpara.bath],'Motion JPEG AVI');
    outputVideo.FrameRate=2;
    open(outputVideo)
    
    for t=2000:2200
    figure(1)
    subplot(211)
    curlF_temp = curlF(:,:,t);
    pcolorcen(y_grid,x_grid,curlF_temp');
    col=colorbar;
    cmocean('balance');
    %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
    caxis([-1.0,1.0])
    col.Label.Interpreter = 'latex';
    col.Label.String = '$\nabla \times F_{br} \, (s^{-2})$';
    col.Label.FontSize = 22;
    col.Label.FontWeight = 'bold';
    ylabel('x (m)','Interpreter','latex')
    xlabel('y (m)','Interpreter','latex')
    title([runpara.wave,runpara.bath])
    width= 25;
    height = 10;
    set(gcf,'Units','inches','Position',[0,0,width,height]);
    set(gcf,'visible','off');
    ylim([1,dim(1)])
    xlim([0,350])
    set(gca,'YDir','reverse')
    niceplot(22)
    
    
    subplot(212)
    vort_temp = vort(:,:,t);
    pcolorcen(y_grid,x_grid,vort_temp');
    col=colorbar;
    cmocean('balance');
    %caxis([-abs(max(abs(cxt_atx(:)),[],'all')),abs(max(abs(cxt_atx(:)),[],'all'))])
    caxis([-0.25,0.25])
    col.Label.Interpreter = 'latex';
    col.Label.String = '$\nabla \times u \, (s^{-1})$';
    col.Label.FontSize = 22;
    col.Label.FontWeight = 'bold';
    ylabel('x (m)','Interpreter','latex')
    xlabel('y (m)','Interpreter','latex')
    title([runpara.wave,runpara.bath])
    width= 22;
    height = 18;
    set(gcf,'Units','inches','Position',[0,0,width,height]);
    set(gcf,'visible','off');
    %ylim([dim(1)-100,dim(1)])
    ylim([1,dim(1)])
    xlim([0,350])
    set(gca,'YDir','reverse')
    niceplot(22)
    
    
    % scale=0.1;
    % pos = get(gca,'Position');
    % pos(1) = pos(2)+scale*pos(4)-0.02;
    % pos(3) = (1-scale)*pos(4);
    % set(gca,'Position',pos)
    
    frame = getframe(gcf);
    writeVideo(outputVideo,frame)
    close(gcf)
    end 
    close(outputVideo)
    runind
    toc
%end 