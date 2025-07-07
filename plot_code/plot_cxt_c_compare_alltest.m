% Bingchen Liu March 11, 2025
% This code compare the linear regression of dx location of max CXT to
% c=\sqrt(gh) prediction to validate the use of choice of max CXT 

clear
load('/data1/bliu/data/cxt_alongct_nointerp_max_dxwidth') %'get_cxt_alongct_nointerp_max_dxwidth'
load('/data1/bliu/data/runnum_72run')



markersz= 50;


%% c/sqrt(gh) VS h/hb
factor = 1;
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

c_nond_diffxloc=zeros(72,5);
for xloc=1:5
for i =1:24
    c_nond_diffxloc(3*i-2:3*i,xloc)= [c_phase_5loc.slp2(i).c_nond(xloc),c_phase_5loc.slp3(i).c_nond(xloc),c_phase_5loc.slp4(i).c_nond(xloc)];

end 
end 

c_stats.mean  = mean(c_nond_diffxloc,1,'omitmissing');
c_stats.mean_all = mean(c_nond_diffxloc(:),1,'omitmissing');
c_stats.std  = std(c_nond_diffxloc,1,'omitmissing');


c_stats_xloc=  [0.25,0.37,0.5,0.62,0.75]; %estimated h/hb for all 5 sz loc
c_stats_xloc=  flip(c_stats_xloc);
figure()
%subplot(121)
for xloc =x_num:-1:1
    for N =1:24
        scatter((c_phase_5loc.slp2(N).h_hb_nond(xloc)).^factor,c_phase_5loc.slp2(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter((c_phase_5loc.slp3(N).h_hb_nond(xloc)).^factor,c_phase_5loc.slp3(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter((c_phase_5loc.slp4(N).h_hb_nond(xloc)).^factor,c_phase_5loc.slp4(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
    end
end 

errorbar(c_stats_xloc,c_stats.mean,c_stats.std)
%plot_oneone
hold off
xlabel('$h/h_b$','Interpreter','latex')
ylabel('$c/ \sqrt{g h}$','Interpreter','latex')
niceplot(12)


%%
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

figure()
%subplot(121)
for xloc =x_num:-1:1
    for N =1:24
        scatter(c_phase_5loc.slp2(N).c_modelh(xloc),c_phase_5loc.slp2(N).c_fit(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp3(N).c_modelh(xloc),c_phase_5loc.slp3(N).c_fit(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp4(N).c_modelh(xloc),c_phase_5loc.slp4(N).c_fit(xloc),markersz,col(xloc,:),'filled')
        hold on 
    end
end 
xlim([0,7])
ylim([0,7])
plot_oneone
hold off
xlabel('$c = \sqrt{g h} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
ylabel('$\hat{c} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
niceplot(12)

% subplot(122)
% for xloc =x_num:-1:1
%     for N =1:24
%         scatter(c_phase_5loc.slp2(N).c_modelhHs(xloc),c_phase_5loc.slp2(N).c_fit(xloc),markersz,col(xloc,:),'filled')
%         hold on 
%         scatter(c_phase_5loc.slp3(N).c_modelhHs(xloc),c_phase_5loc.slp3(N).c_fit(xloc),markersz,col(xloc,:),'filled')
%         hold on 
%         scatter(c_phase_5loc.slp4(N).c_modelhHs(xloc),c_phase_5loc.slp4(N).c_fit(xloc),markersz,col(xloc,:),'filled')
%         hold on 
%     end
% end 
% xlim([0,7])
% ylim([0,7])
% plot_oneone
% hold off
% xlabel('$c = \sqrt{g (h+H_{s})} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% ylabel('$\hat{c} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% niceplot(12)
% 
% clear col
% colormap(cmocean('thermal',x_num))
% cbar=colorbar('south');
% clim([-0.75-0.125/2,-0.25+0.125/2]);
% cbar.Ticks =linspace(-0.75,-0.25,5);
% cbar.Label.FontSize = 14;
% cbar.Label.FontWeight = 'bold';
% cbar.Label.Interpreter = 'latex';
% cbar.Label.String = '$x / L_{\mathrm{sz}}$';
% %cbar.Position = [0.3, 0.12,0.5,0.02];
%%
x_num = length(c_phase_5loc.slp2(1).c_fit);
col = cmocean('thermal',x_num);

figure()
subplot(121)
for xloc =x_num:-1:1
    for N =1:24
        scatter(c_phase_5loc.slp2(N).Hs_h(xloc),c_phase_5loc.slp2(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp3(N).Hs_h(xloc),c_phase_5loc.slp3(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp4(N).Hs_h(xloc),c_phase_5loc.slp4(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
    end
end 
%xlim([0,7])
%ylim([0,7])
%plot_oneone
hold off
xlabel('$\frac{H_{s}}{h} $','Interpreter','latex')
ylabel('$\hat{c}/c $','Interpreter','latex')
niceplot(12)

subplot(122)
for xloc =x_num:-1:1
    for N =1:24
        scatter(c_phase_5loc.slp2(N).Hs_h_sqrt(xloc),c_phase_5loc.slp2(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp3(N).Hs_h_sqrt(xloc),c_phase_5loc.slp3(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
        scatter(c_phase_5loc.slp4(N).Hs_h_sqrt(xloc),c_phase_5loc.slp4(N).c_nond(xloc),markersz,col(xloc,:),'filled')
        hold on 
    end
end 
%xlim([0,7])
%ylim([0,7])
%plot_oneone
hold off
xlabel('$ \sqrt{1+\frac{H_{s}}{h}} $','Interpreter','latex')
ylabel('$\hat{c}/c $','Interpreter','latex')
niceplot(12)

clear col
colormap(cmocean('thermal',x_num))
cbar=colorbar('south');
clim([-0.75-0.125/2,-0.25+0.125/2]);
cbar.Ticks =linspace(-0.75,-0.25,5);
cbar.Label.FontSize = 14;
cbar.Label.FontWeight = 'bold';
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$x / L_{\mathrm{sz}}$';


%% test nond_c VS directional spread
figure()
for xloc =x_num:-1:1
    for N =1:24
        scatter(c_phase_5loc.slp2(N).dirspr_b(xloc),c_phase_5loc.slp2(N).c_nond(xloc),markersz,"red",'filled')
        hold on 
        scatter(c_phase_5loc.slp3(N).dirspr_b(xloc),c_phase_5loc.slp3(N).c_nond(xloc),markersz,"green",'filled')
        hold on 
        scatter(c_phase_5loc.slp4(N).dirspr_b(xloc),c_phase_5loc.slp4(N).c_nond(xloc),markersz,'blue','filled')
        hold on 
    end
end 
%xlim([0,7])
%ylim([0,7])
%plot_oneone
hold off
xlabel('$ \sigma_{\theta b} $','Interpreter','latex')
ylabel('$\hat{c}/c $','Interpreter','latex')
niceplot(12)

%%
% figure()
% subplot(121)
% for i =1:72
%     runind = goodrunnum(i);
%     scatter(c_modelh{runind},c_fit{runind})
%     hold on 
% end 
% xlim([0,10])
% ylim([0,10])
% xlabel('$c = \sqrt{gh} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% ylabel('$\hat{c} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% plot_oneone
% hold off 
% niceplot(12)
% 
% subplot(122)
% for i =1:72
%     runind = goodrunnum(i);
%     scatter(c_modelhHs{runind},c_fit{runind})
%     hold on 
% end  
% xlim([0,10])
% ylim([0,10])
% plot_oneone
% hold off
% xlabel('$c = \sqrt{g (h+H_{s})} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% ylabel('$\hat{c} \, ( \mathrm{m} \,\mathrm{s}^{-1})$','Interpreter','latex')
% niceplot(12)
