% Bingchen Liu June 27, 2024
% This code compare curlFbr (true model data output) and G0
% (inetegration from Weibull fit) 

% Note this is in function format, and is included in master_plot_withloop
% to get plot for each/all run

% NOTE: it is still in progress July 5, 2024

function plot_runnum_curlFbr_Go_compare(ind)

    close all
    addpath('/data1/nkumar/RIPX/M_Files')
    addpath('/data1/bliu')
    
    load('/data1/bliu/data/RIPX_bath_guide.mat')
    load('/data1/bliu/data/SS_raw.mat')
    load('/data1/bliu/data/SS_A_ds_qced_3loc')

    
    %% grab data
    runpar= get_runpara(ind);

    
    [p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02); %p2 is runnumber; v2=[Hs3 fp3 m3 spread3];
    [p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
    [p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);


    G0_wb = zeros(120,3);
    
    index_slpgroup = [p2;p3;p4];
    G0_wb_temp= [A.wfit2.intw;A.wfit3.intw;A.wfit4.intw];
    
    
    for  i = 1:length(G0_wb_temp)
        G0_wb(index_slpgroup(i),:) = G0_wb_temp(i,:);
    end
    
    
    
    x= S(ind).X;
    x2 = S(ind).X2;   
    
    ind_3loc_1 = get_3locs(S(ind),1);
    ind_3loc_2 = get_3locs(S(ind),2);    
    
    curlFstd = S(ind).curlF_std;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
    %xrange=length(x2)-10:length(x2); %plot only section of the scaling to prevent sqrt(gh) blowing up 
    figure()
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20)

    plot(x,curlFstd,'Linewidth',3)
    %hold on 
    %scatter(x(ind_3loc_1),G0_wb(ind,:),80,"filled")
    %hold off 
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$std(\nabla \times F_{br}) \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    niceplot_nobold_nomintick(22)
    grid on   
    
    figure()
    sgtitle([runpar.wave,runpar.bath],'Fontsize',20) 
    scatter(x(ind_3loc_1),G0_wb(ind,:),80,"filled")
    xlabel('Cross-shore[x(m)]','Fontsize',16)
    ylabel('$std(\nabla \times F_{br}) \, (s^{-2})$','Interpreter','latex')
    xlim([min(x2),max(x2)])
    ylim([min(G0_wb(ind,:)),max(G0_wb(ind,:))])
    niceplot_nobold_nomintick(22)
    grid on   
    
end 