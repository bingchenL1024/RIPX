% Bingchen Liu June 25, 2025
% This code compare different scaling for ky0


clear
load('/data1/bliu/data/plotready_G0_nond_10loc_2025')
load('/data1/bliu/data/G0_nond_10loc_allsk_2025')


%%%%%%%%%%%%%%%%%%%%%%%%%
loc_num =length(G0_nond.slp2(1,:));

sk_ky0.Sk2_h=1./(SS.h2);
sk_ky0.Sk3_h=1./(SS.h3);
sk_ky0.Sk4_h=1./(SS.h4);


sk_ky0.Sk2_h_beta=(1./SS.h2)./SS.beta2;
sk_ky0.Sk3_h_beta=(1./SS.h3)./SS.beta3;
sk_ky0.Sk4_h_beta=(1./SS.h4)./SS.beta4;

sk_ky0.Sk2_h_Ir=(1./SS.h2)./SS.Irb2;
sk_ky0.Sk3_h_Ir=(1./SS.h3)./SS.Irb3;
sk_ky0.Sk4_h_Ir=(1./SS.h4)./SS.Irb4;

ky0_nond.slp2 = A.wkym2./sk_ky0.Sk2_h;
ky0_nond.slp3 = A.wkym3./sk_ky0.Sk3_h;
ky0_nond.slp4 = A.wkym4./sk_ky0.Sk4_h;

ky0_nond_beta.slp2 = A.wkym2./sk_ky0.Sk2_h_beta;
ky0_nond_beta.slp3 = A.wkym3./sk_ky0.Sk3_h_beta;
ky0_nond_beta.slp4 = A.wkym4./sk_ky0.Sk4_h_beta;


ky0_nond_Ir.slp2 = A.wkym2./sk_ky0.Sk2_h_Ir;
ky0_nond_Ir.slp3 = A.wkym3./sk_ky0.Sk3_h_Ir;
ky0_nond_Ir.slp4 = A.wkym4./sk_ky0.Sk4_h_Ir;



%% calculate rsquare 
G0=load('/data1/bliu/data/plotready_G0_nond_10loc_2025.mat');

loc_num=length(G0.G0_nond.slp2(1,:));
x_num = length(G0_nond.slp2(1,:));

ds_10loc.slp2 =repmat(G0.ds.sigtb2,[loc_num,1])';
ds_10loc.slp3 =repmat(G0.ds.sigtb3,[loc_num,1])';
ds_10loc.slp4 =repmat(G0.ds.sigtb4,[loc_num,1])';

rsq = get_rsq(ky0_nond,ds_10loc)
rsq_beta = get_rsq(ky0_nond_beta,ds_10loc)
rsq_Ir = get_rsq(ky0_nond_Ir,ds_10loc)


%% plot 
markersz=35;


col = cmocean('thermal',x_num);


figure()
subplot(131)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0_nond.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),ky0_nond.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),ky0_nond.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$$ \hat{k}_{y0} \, h$$','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)

subplot(132)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0_nond_beta.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),ky0_nond_beta.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),ky0_nond_beta.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$$ \hat{k}_{y0} \, h \, \beta$$','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)

subplot(133)
for xind=1:x_num
scatter(ds_10loc.slp2(:,xind),ky0_nond_Ir.slp2(:,xind),markersz,col(xind,:),'filled','Marker','o','MarkerEdgeColor','k')
scatter(ds_10loc.slp3(:,xind),ky0_nond_Ir.slp3(:,xind),markersz,col(xind,:),'filled','Marker','square','MarkerEdgeColor','k')
scatter(ds_10loc.slp4(:,xind),ky0_nond_Ir.slp4(:,xind),markersz,col(xind,:),'filled','Marker','^','MarkerEdgeColor','k')
hold on 
end 
hold off
ylabel('$$ \hat{k}_{y0} \, h \, Ir$$','interpreter','latex');
xlabel('$\sigma_{\theta b}$~(deg)','interpreter','latex');
niceplot(26)


%% functions 


function rsq = get_rsq(G0_nond,ds)

G0_nond_all = [G0_nond.slp2(:);G0_nond.slp3(:);G0_nond.slp4(:);];
ds_all = [ds.slp2(:),ds.slp3(:),ds.slp4(:)];
ftt = strcat('A*x+B');
ft = fittype( sprintf('%s',ftt));
opts = fitoptions( ft );
opts.Display = 'Off';
%opts.Lower = [0 -0.5];
opts.StartPoint = [0.5 2]; % beginning parameters - amp, mu, std.
%opts.Upper = [Inf Inf]; %phase shift should be within 2pi

[f,gof]=fit(ds_all(~isnan(G0_nond_all)),G0_nond_all(~isnan(G0_nond_all)),ft, opts);
G0_pred = f.A*ds_all(~isnan(G0_nond_all))+f.B;

rsq= falk_skill(G0_pred,G0_nond_all(~isnan(G0_nond_all)));
%rsq=gof.rsquare;
end 