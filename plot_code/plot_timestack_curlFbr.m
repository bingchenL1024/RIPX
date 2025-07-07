% Bingchen Liu Jan 6, 2025 
% This code plot time stack of curlFbr at given y location and for a given
% run

clear all 
%%
runind=44;

clearvars -except runind
fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',runind);
load(fname,'curlF')
dim_curlFbr= size(curlF);


yind = 200;
tlim = 200:350;
fontsz = 26;

curlF_tstack = squeeze(curlF(:,yind,tlim));
save('/data1/bliu/data/tstack_curlF_run44', "curlF_tstack")
%%
x = -(dim_curlFbr(1)-1):0;
t = 0:length(tlim)-1; %in seconds
[t_grid,x_grid] = meshgrid(t,x);

pcolorcen(t_grid,x_grid,tstack);
col=colorbar;
cmocean('balance');
caxis([-0.5,0.5])
set(gca,'YDir','reverse') %modified

col.Label.Interpreter = 'latex';
col.Label.String = '$\nabla \times \mathbf{F}_{\mathrm{br}} \, (s^{-2})$';
col.Label.FontSize = fontsz;

ylabel('$x$ (m)','Interpreter','latex','FontSize',fontsz)
xlabel('$t$ (m)','Interpreter','latex','FontSize',fontsz)

ylim([-150,0]) %modified

title(['run',num2str(runind),', y =',num2str(yind), ', t = ' , num2str(tlim(1)),':',num2str(tlim(end))],'FontSize',fontsz);



















