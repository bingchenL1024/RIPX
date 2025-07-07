% Bingchen Liu July 5, 2024
% This code plots raw spectra (no Weibull fit) from model output raw
% spectra
% Originally designed to debug Weibull ky0 misfit at locations that are
% very close to the shoreline -- > 


%ind_run: index of the runs that you want to work with (basically runnumber)
%ind_plot: gives locations where you want the spectra plot. It is such that
%X location(-XXX m): C(ind_plot(1)) --> give you X location  

function plot_runnum_spectra_at_diffloc(ind_run)


load('/data1/bliu/data/SS_raw.mat')

runpar= get_runpara(ind_run);

SS4 = S(ind_run);
h = SS4.h;

om = SS4.fp_T.*(2.*pi);
kw = get_wavenum(om,h); % use omega at breaking pt 
kw_pm = (kw./(2*pi));

kh = h.* kw_pm;

P4 = load_RIPX_curlF_Sky(ind_run);
Sky4 = P4.Sckm; %data spectra dim: crosshore*spectra dim
kyy = P4.fck/(2*pi); %data
[C,IA,IB] = intersect(SS4.X,SS4.X2); %29 is the cross-shore dimension 
Sky4 = Sky4(IA,:); % now Sky4 has the same dimension as C


if length(C)> 30;  % Note: 1. C(ind_plot) --> gives X location 2. it only plots the loc that is very close to shore
    ind_plot = length(C)-16:4:length(C);
elseif length(C) > 20;
    ind_plot = length(C)-12:3:length(C); 
elseif length(C)>13;
    ind_plot = length(C)-8:2:length(C); 
else length(C)<=13;;
    ind_plot = length(C)-4:1:length(C); 
end 

% 
% 
% if length(C)<9
%     ind_plot = length(C)-4:1:length(C); 
% elseif length(C)<17%set cross shore locations that will be plotted (based on the domain)
%     ind_plot = length(C)-12:3:length(C); 
% else 
%     ind_plot = length(C)-16:4:length(C); 
% end 

%%%%%%%%%%%%%%%%%%%%%%% doing Weibull fit 
%for xx = 1:length(IA)
for xx = 1:length(ind_plot)

%%%%%%%%%%%%%%%%%% interpolate original raw data spectra
[S4] = interpSky(Sky4(ind_plot(xx),:),kyy,91);
%iS4s(xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));
%df = diff(S4.lky);
lky(xx,:) = S4.lky;
% Skylo(xx,:) = S4.Skyl;
% Skyls(xx,:) = S4.sSkyl;


%%%%%%%%%%%%%%%%% Do the Weibull fit:
[FP4] = rayleigh_fitAS2(S4.lky,S4.sSkyl,1.2);
% lskil4(xx) = falk_skill(log10(FP4.line'),log10(S4.sSkyl));
% wskil4(xx) = wilmot_skill(log10(FP4.line'),log10(S4.sSkyl));
%intw(xx) = nansum(FP4.line.*[0 df]');
wSkyl(xx,:) = FP4.line;

a = find(S4.sSkyl==max(S4.sSkyl));
Skym4(xx) = S4.sSkyl(a);
kym4(xx) = S4.lky(a);
a = find(FP4.line==max(FP4.line));
wSkym4(xx) = FP4.line(a);
wkym4(xx) = S4.lky(a);
end







kyy_plot = repmat(kyy',length(ind_plot),1);
colormap = parula(length(ind_plot));
figure()
sgtitle([runpar.wave,runpar.bath],'Fontsize',20)
for i = 1:length(ind_plot)
loglog(kyy_plot(i,:)',Sky4(ind_plot(i),:)','linewidth',3,'Color',colormap(i,:),'DisplayName',['h = ', num2str(h(IA(ind_plot(i)))),' m ', ' kh = ',num2str(kh(IA(ind_plot(i))),sprintf('%%.%df',3))]) % pre-interpolation
hold on 
loglog(lky(i,:)',wSkyl(i,:),'--','color',[.5 .5 .5],'linewidth',2,'Color',colormap(i,:),'DisplayName','Weibull Fit')
end 
legend show
legend('Location','southwest')
ylabel('$S_{\nabla \times \mathrm{Fbr}}(k_y) \ \ (\mathrm{s^{-4}/cpm})$','interpreter','latex','fontsize',16);
xlabel('$k_{y} \ (\mathrm{cpm})$','interpreter','latex','fontsize',18);
xlim([0,0.2])
ylim([min(Sky4(ind_plot,:),[],'all'),max(Sky4(ind_plot,:),[],'all')])
niceplot_nobold_nomintick(22);
hold off

width= 18;
height = 15;
set(gcf,'Units','inches','Position',[0,0,width,height])
set(gcf,'visible','off') 

saveas(gcf,['/data1/bliu/figures/RIPX_allrun/Spectra/','run_',num2str(ind_run),'.png'])

end 