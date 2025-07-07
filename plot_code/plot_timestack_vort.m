% Bingchen Liu Jan 22, 2025 
% This code plot time stack of curlFbr at given y location and for a given
% run

clear all 
%%
runind=44;

clearvars -except runind
fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_vort/RIPX_vort_%04d',runind);

load(fname,'vort')
%vort = permute(vort,[2 3 1]);
dim_vort= size(vort);


yind = 200;
tlim = 200:350;
fontsz = 26;

vort_tstack = squeeze(vort(tlim,:,yind))';
save('/data1/bliu/data/tstack_vort_run44', "vort_tstack")
%%
x = -(dim_vort(2)-1):0;
t = 0:length(tlim)-1; %in seconds
[t_grid,x_grid] = meshgrid(t,x);

pcolorcen(t_grid,x_grid,tstack);
shading interp
col=colorbar;
cmocean('balance');
caxis([-0.5,0.5])
set(gca,'YDir','reverse') %modified

col.Label.Interpreter = 'latex';
col.Label.String = '$\omega \, (s^{-1})$';
col.Label.FontSize = fontsz;

ylabel('$x$ (m)','Interpreter','latex','FontSize',fontsz)
xlabel('$t$ (m)','Interpreter','latex','FontSize',fontsz)

ylim([-150,0]) %modified

title(['run',num2str(runind),', y =',num2str(yind), ', t = ' , num2str(tlim(1)),':',num2str(tlim(end))],'FontSize',fontsz);


