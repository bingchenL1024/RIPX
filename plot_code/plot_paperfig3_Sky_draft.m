
% Bingchen Liu Nov 22, 2024
% This code is used to plot along-shore wavenumber spectra 
% Correspond to Fig 5 without the Weibull fit,
% NOTE: this code uses specific run, for general spec plot, check 

clear
load('/data1/bliu/data/Sky_nond')
load('/data1/bliu/data/Sky_WBfit_qced') 
load('/data1/bliu/data/plotready_ky0_nond_3loc.mat')


figure()
% subplot(121)
for N= 1:24
loglog(squeeze(ky2_cutoff(N,:,:))',squeeze(Sky2_cutoff_runmean(N,:,:))')
hold on 
loglog(squeeze(ky3_cutoff(N,:,:))',squeeze(Sky3_cutoff_runmean(N,:,:))')
hold on 
loglog(squeeze(ky4_cutoff(N,:,:))',squeeze(Sky4_cutoff_runmean(N,:,:))')

end 
hold off 
ylabel('$S_{\nabla \times F_{br}}(s^{-4}/cpm)$','Interpreter','latex')
xlabel('$k_{y} (cpm)$','Interpreter','latex')
niceplot_nobold_nomintick(22)
% 
% subplot(132)
% for N = 1:24
% loglog(squeeze(ky2_cutoff(N,:,:))',squeeze(Sky2_nond_maxsky(N,:,:))')
% hold on 
% end 
% ylabel('$\hat{S}_{\nabla \times F_{br}}$','Interpreter','latex')
% xlabel('$k_{y} (cpm)$','Interpreter','latex')
% niceplot_nobold_nomintick(22)

% 
% subplot(122)
% for N =1:24
% loglog(squeeze(ky2_nond(N,:,:))',squeeze(Sky2_nond_maxsky(N,:,:))')
% hold on 
% end 
% ylabel('$\hat{S}_{\nabla \times F_{br}}$','Interpreter','latex')
% xlabel('$\hat{k}_{y}$','Interpreter','latex')
% niceplot_nobold_nomintick(22)
%%
figure()

% subplot(131)
% for N= 1:24
% loglog(squeeze(ky2_cutoff(N,1,:))',squeeze(Sky2_nond(N,1,:))','r')
% hold on 
% end 
% for N= 1:24
% loglog(squeeze(ky2_cutoff(N,2,:))',squeeze(Sky2_nond(N,2,:))','g')
% hold on 
% end 
% for N= 1:24
% loglog(squeeze(ky2_cutoff(N,3,:))',squeeze(Sky2_nond(N,3,:))','b')
% hold on 
% end 
% 
% hold off 
% ylabel('$\hat{S}_{\nabla \times F_{br}}$','Interpreter','latex')
% xlabel('$k_{y} (cpm)$','Interpreter','latex')
% niceplot_nobold_nomintick(22)
% 



subplot(131)
for N =1:24
loglog(squeeze(ky2_nond(N,1,:))',squeeze(Sky2_nond_maxsky(N,1,:))','r')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,2,:))',squeeze(Sky2_nond_maxsky(N,2,:))','g')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,3,:))',squeeze(Sky2_nond_maxsky(N,3,:))','b')
hold on 
end
ylabel('$S_{\nabla \times F_{br}}/S_{\mathrm{max}}$','Interpreter','latex')
xlabel('$\hat{k}_{y}$','Interpreter','latex')
xlim([-inf 10])
set(gca,'XTick',[0.01,0.1,1,10])
niceplot_nobold_nomintick(22)
grid on 


subplot(132)
for N =1:24
loglog(squeeze(ky2_nond(N,1,:))',squeeze(Sky2_nond_kym(N,1,:))','r')
hold on
loglog(squeeze(ky3_nond(N,1,:))',squeeze(Sky3_nond_kym(N,1,:))','r')
hold on 
loglog(squeeze(ky4_nond(N,1,:))',squeeze(Sky4_nond_kym(N,1,:))','r')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,2,:))',squeeze(Sky2_nond_kym(N,2,:))','g')
hold on
loglog(squeeze(ky3_nond(N,2,:))',squeeze(Sky3_nond_kym(N,2,:))','g')
hold on 
loglog(squeeze(ky4_nond(N,2,:))',squeeze(Sky4_nond_kym(N,2,:))','g')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,3,:))',squeeze(Sky2_nond_kym(N,3,:))','b')
hold on 
loglog(squeeze(ky3_nond(N,3,:))',squeeze(Sky3_nond_kym(N,3,:))','b')
hold on 
loglog(squeeze(ky4_nond(N,3,:))',squeeze(Sky4_nond_kym(N,3,:))','b')
hold on 
enduntitled5
ylabel('$S_{\nabla \times F_{br}} k_{y0}/ {\int_{0.001}^{0.2} \, S_{\nabla \times F_{br}} \, dk_{y}}$','Interpreter','latex')
xlabel('$\hat{k}_{y}$','Interpreter','latex')
%xlim([-inf 20])
niceplot_nobold_nomintick(22)
end 


for N =1:24
loglog(squeeze(ky2_nond_wb(N,1,:))',squeeze(Sky2_nond_kym_wb(N,1,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky3_nond_wb(N,1,:))',squeeze(Sky3_nond_kym_wb(N,1,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on
loglog(squeeze(ky4_nond_wb(N,1,:))',squeeze(Sky4_nond_kym_wb(N,1,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on
end

for N =1:24
loglog(squeeze(ky2_nond_wb(N,2,:))',squeeze(Sky2_nond_kym_wb(N,2,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky3_nond_wb(N,2,:))',squeeze(Sky3_nond_kym_wb(N,2,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky4_nond_wb(N,2,:))',squeeze(Sky4_nond_kym_wb(N,2,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond_wb(N,3,:))',squeeze(Sky2_nond_kym_wb(N,3,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky3_nond_wb(N,3,:))',squeeze(Sky3_nond_kym_wb(N,3,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky4_nond_wb(N,3,:))',squeeze(Sky4_nond_kym_wb(N,3,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
end
ylabel('$S_{ky} k_{y0}/ {\int_{0.001}^{0.2} \, S_{ky} \, dk_{y}}$','Interpreter','latex')
xlabel('$\hat{k}_{y}$','Interpreter','latex')
%xlim([-inf 20])

%ylabel('$S_{WB} k_{y0}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
%xlabel('$\hat{k}_{y}$','Interpreter','latex')
%xlim([-inf 20])
ylim([10e-3 10e-1])
set(gca,'XTick',[0.01,0.1,1,10])
niceplot_nobold_nomintick(22)
grid on 

subplot(133)
for N =1:24
loglog(squeeze(ky2_nond(N,1,:))',squeeze(Sky2_nond_kym_S0(N,1,:))','r')
hold on 
loglog(squeeze(ky3_nond(N,1,:))',squeeze(Sky3_nond_kym_S0(N,1,:))','r')
hold on 
loglog(squeeze(ky4_nond(N,1,:))',squeeze(Sky4_nond_kym_S0(N,1,:))','r')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,2,:))',squeeze(Sky2_nond_kym_S0(N,2,:))','g')
hold on
loglog(squeeze(ky3_nond(N,2,:))',squeeze(Sky3_nond_kym_S0(N,2,:))','g')
hold on 
loglog(squeeze(ky4_nond(N,2,:))',squeeze(Sky4_nond_kym_S0(N,2,:))','g')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond(N,3,:))',squeeze(Sky2_nond_kym_S0(N,3,:))','b')
hold on 
loglog(squeeze(ky3_nond(N,3,:))',squeeze(Sky3_nond_kym_S0(N,3,:))','b')
hold on 
loglog(squeeze(ky4_nond(N,3,:))',squeeze(Sky4_nond_kym_S0(N,3,:))','b')
hold on 
end

for N =1:24
loglog(squeeze(ky2_nond_wb(N,:,:))',squeeze(Sky2_nond_kym_wb_S0(N,:,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on
loglog(squeeze(ky3_nond_wb(N,:,:))',squeeze(Sky3_nond_kym_wb_S0(N,:,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
loglog(squeeze(ky4_nond_wb(N,:,:))',squeeze(Sky4_nond_kym_wb_S0(N,:,:))','LineWidth',5,'Color',[.5 .5 .5])
hold on 
end
ylabel('$S_{\nabla \times F_{br}} k_{y0\_WB}/ {\int_{0}^{\infty} \, S_{WB} \, dk_{y}}$','Interpreter','latex')
xlabel('$\hat{k}_{y}$','Interpreter','latex')
%xlim([-inf 10])
set(gca,'XTick',[0.01,0.1,1,10])
niceplot_nobold_nomintick(22)
grid on 
axis([0.02 10 0.01 1])

% %% test&debug
% 
% for N = 1:24
% loglog(squeeze(ky2_nond_wb(N,:,:))',squeeze(Sky2_nond_kym_wb_S0(N,:,:))','LineWidth',5,'Color',[.5 .5 .5])
% hold on
% end 
% hold off
% 
% 
% test_ky = squeeze(ky2_nond_wb(3,:,:));
% test_sky = squeeze(Sky2_nond_kym_wb_S0(3,:,:));
% 
% test_ky1 = squeeze(ky2_nond_wb(2,:,:));
% test_sky1 = squeeze(Sky2_nond_kym_wb_S0(2,:,:));
% 
% figure()
% for i = 1:3
% loglog(test_ky(i,:),test_sky(i,:));hold on;loglog(test_ky1(i,:),test_sky1(i,:));hold off
% end 