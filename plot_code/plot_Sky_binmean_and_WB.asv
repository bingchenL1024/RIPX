% Bingchen Liu Jan 15, 2025
% This code plot the nond bin meaned Sky along with collapsed Weibll fit 
% first figure only shows one plot of given exp; second figure shows 3 exp
% for comparison (currently disabled)

clearvars -except expo

expo_name_temp = num2str(expo);
expo_name  = [expo_name_temp(1),'p',expo_name_temp(3:end)];

load(['/data1/bliu/data/Sky_binmean_',expo_name])



figure()
errorbar(ky_bincenter,Sky_binmean,error_binmean,'s','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','auto')
set(gca,'XScale','log','YScale','log')

hold on 
plot(ky_nond,Sky_nond_analyt,'LineWidth',3,'Color','r','LineStyle','--')

% for N =1:24
% h2=loglog(squeeze(ky2_nond_wb(N,:,:))',squeeze(Sky2_nond_kym_wb_S0(N,:,:))','LineWidth',3,'Color','r','LineStyle','--');
% hold on
% loglog(squeeze(ky3_nond_wb(N,:,:))',squeeze(Sky3_nond_kym_wb_S0(N,:,:))','LineWidth',3,'Color','r','LineStyle','--')
% hold on 
% loglog(squeeze(ky4_nond_wb(N,:,:))',squeeze(Sky4_nond_kym_wb_S0(N,:,:))','LineWidth',3,'Color','r','LineStyle','--')
% hold on 
% end

% plot bin edge 
% for i = 1:length(edge_bin)
%     xline(10.^(edge_bin(i)))
% end 

ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
xlabel('$k_{y}/k_{y0}$','Interpreter','latex')
title(['Expo = ',expo_name,', r^2=',num2str(R_square)])
%xlim([-inf 10])
%title(['Expo = ',expo_name])
set(gca,'XTick',[0.001,0.01,0.1,1,10])
axis([0.03 10 0.01 1])
niceplot_nobold(23)
grid on 

%% plot diff expo on the same plot for comparison 
% clear
% expo_all = 1.225:0.025:1.275;
% 
% for i=1:length(expo_all)
%     expo = expo_all(i);
%     expo_name_temp = num2str(expo);
%     expo_name  = [expo_name_temp(1),'p',expo_name_temp(3:end)];
%     diffexp = ['exp_',num2str(i)];
%     A.(diffexp)=load(['/data1/bliu/data/Sky_binmean_',expo_name]);
%     A.(diffexp).expo = expo;
% end 
% 
% 
% figure()
% subplot(131)
% errorbar(A.exp_1.ky_bincenter,A.exp_1.Sky_binmean,A.exp_1.error_binmean,'s','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','auto')
% set(gca,'XScale','log','YScale','log')
% 
% hold on 
% plot(A.exp_1.ky_nond,A.exp_1.Sky_nond_analyt,'LineWidth',3,'Color','r','LineStyle','--')
% 
% ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
% xlabel('$k_{y}/k_{y0}$','Interpreter','latex')
% %xlim([-inf 10])
% title(['Expo = ',num2str(A.exp_1.expo),', r^2=',num2str(A.exp_1.R_square),', RMSE=',num2str(A.exp_1.rmse_bin)],'FontWeight','normal')
% set(gca,'XTick',[0.001,0.01,0.1,1,10])
% axis([0.03 10 0.01 1])
% niceplot_nobold(16)
% grid on 
% hold off
% 
% 
% subplot(132)
% errorbar(A.exp_2.ky_bincenter,A.exp_2.Sky_binmean,A.exp_2.error_binmean,'s','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','auto')
% set(gca,'XScale','log','YScale','log')
% hold on 
% plot(A.exp_2.ky_nond,A.exp_2.Sky_nond_analyt,'LineWidth',3,'Color','r','LineStyle','--')
% 
% ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
% xlabel('$k_{y}/k_{y0}$','Interpreter','latex')
% %xlim([-inf 10])
% title(['Expo = ',num2str(A.exp_2.expo),', r^2=',num2str(A.exp_2.R_square),', RMSE=',num2str(A.exp_2.rmse_bin)],'FontWeight','normal')
% set(gca,'XTick',[0.001,0.01,0.1,1,10])
% axis([0.03 10 0.01 1])
% niceplot_nobold(16)
% grid on 
% hold off 
% 
% subplot(133)
% errorbar(A.exp_3.ky_bincenter,A.exp_3.Sky_binmean,A.exp_3.error_binmean,'s','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','auto')
% set(gca,'XScale','log','YScale','log')
% hold on 
% plot(A.exp_3.ky_nond,A.exp_3.Sky_nond_analyt,'LineWidth',3,'Color','r','LineStyle','--')
% 
% ylabel('$S_{\nabla \times \mathbf{F}_{\mathrm{br}}} \, k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
% xlabel('$k_{y}/k_{y0}$','Interpreter','latex')
% %xlim([-inf 10])
% title(['Expo = ',num2str(A.exp_3.expo),', r^2=',num2str(A.exp_3.R_square),', RMSE=',num2str(A.exp_3.rmse_bin)],'FontWeight','normal')
% set(gca,'XTick',[0.001,0.01,0.1,1,10])
% axis([0.03 10 0.01 1])
% niceplot_nobold(16)
% grid on 
% hold off 
