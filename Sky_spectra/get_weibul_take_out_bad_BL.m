% Skill of the fit. Skill done in Log space:
addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))
clear
load /data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat % now S(1--120) struct is in memory
load('/data1/bliu/data/Sky_withWBfit')

%%
expo=1.2;
expo_name_temp = num2str(expo);
expo_name  = [expo_name_temp(1),'p',expo_name_temp(3)];

A = load(['/data1/bliu/data/Sky_WBfit_expo_',expo_name]);
SS = load('Syscale');
K = load('locky');

[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

% Get rid of skills that are less than 0
bla2 = find(A.lskil2(:,1)<=0.0);
bla3 = find(A.lskil3(:,1)<=0.0);
bla4 = find(A.lskil4(:,1)<=0.0);

bla22 = find(A.lskil2(:,2)<=0.0);
bla32 = find(A.lskil3(:,2)<=0.0);
bla42 = find(A.lskil4(:,2)<=0.0);

% Also, have a run in the 0.02 slope that has Weibul ky too large:
[I2,J2] = ind2sub([24 3],50);


%%
% indexs of those we toss:
indf.i2 = unique([bla2;bla22;I2]);
indf.i3 = unique([bla3;bla32]);
indf.i4 = unique([bla4;bla42]);
parsf.dirsp2 = v2(indf.i2,4);
parsf.dirsp3 = v3(indf.i3,4);
parsf.dirsp4 = v4(indf.i4,4);
parsf.hs2 = v2(indf.i2,1);
parsf.hs3 = v3(indf.i3,1);
parsf.hs4 = v4(indf.i4,1);
parsf.tp2 = v2(indf.i2,2);
parsf.tp3 = v3(indf.i3,2);
parsf.tp4 = v4(indf.i4,2);

% Resave everything:
K.h2(indf.i2,:) = NaN;K.h3(indf.i3,:) = NaN;K.h4(indf.i4,:) = NaN;
K.kl2(indf.i2,:) = NaN;K.kl3(indf.i3,:) = NaN;K.kl4(indf.i4,:) = NaN;
K.kl2pm(indf.i2,:) = NaN;K.kl3pm(indf.i3,:) = NaN;K.kl4pm(indf.i4,:) = NaN;

A.iS2tot(indf.i2,:) = NaN;A.iS3tot(indf.i3,:) = NaN;A.iS4tot(indf.i4,:) = NaN;
A.iS2s(indf.i2,:) = NaN;A.iS3s(indf.i3,:) = NaN;A.iS4s(indf.i4,:) = NaN;
A.lskil2(indf.i2,:) = NaN;A.lskil3(indf.i3,:) = NaN;A.lskil4(indf.i4,:) = NaN;
A.wskil2(indf.i2,:) = NaN;A.wskil3(indf.i3,:) = NaN;A.wskil4(indf.i4,:) = NaN;

A.Skym2(indf.i2,:) = NaN;A.Skym3(indf.i3,:) = NaN;A.Skym4(indf.i4,:) = NaN;
A.wSkym2(indf.i2,:) = NaN;A.wSkym3(indf.i3,:) = NaN;A.wSkym4(indf.i4,:) = NaN;
A.kym2(indf.i2,:) = NaN;A.kym3(indf.i3,:) = NaN;A.kym4(indf.i4,:) = NaN;
A.wkym2(indf.i2,:) = NaN;A.wkym3(indf.i3,:) = NaN;A.wkym4(indf.i4,:) = NaN;
A.wfit2.a(indf.i2,:) = NaN;A.wfit3.a(indf.i3,:) = NaN;A.wfit4.a(indf.i4,:) = NaN;
A.wfit2.b(indf.i2,:) = NaN;A.wfit3.b(indf.i3,:) = NaN;A.wfit4.b(indf.i4,:) = NaN;
A.wfit2.intw(indf.i2,:) = NaN;A.wfit3.intw(indf.i3,:) = NaN;A.wfit4.intw(indf.i4,:) = NaN;
A.wfit2.So(indf.i2,:) = NaN;A.wfit3.So(indf.i3,:) = NaN;A.wfit4.So(indf.i4,:) = NaN;

SS.h2(indf.i2,:) = NaN; SS.h3(indf.i3,:) = NaN;SS.h4(indf.i4,:) = NaN;
SS.kl2(indf.i2,:) = NaN;  SS.kl3(indf.i3,:) = NaN;SS.kl4(indf.i4,:) = NaN;
SS.Dw2(indf.i2,:) = NaN; SS.Dw3(indf.i3,:) = NaN;SS.Dw4(indf.i4,:) = NaN;
SS.kl2pm(indf.i2,:) = NaN; SS.kl3pm(indf.i3,:) = NaN;SS.kl4pm(indf.i4,:) = NaN;


Sky2_cutoff_runmean(indf.i2,:,:) = NaN; 
Sky3_cutoff_runmean(indf.i3,:,:) = NaN; 
Sky4_cutoff_runmean(indf.i4,:,:) = NaN; 

magspec2(indf.i2,:)  = NaN;
magspec3(indf.i3,:)  = NaN;
magspec4(indf.i4,:)  = NaN;

Skym2(indf.i2,:) = NaN;
Skym3(indf.i3,:) = NaN;
Skym4(indf.i4,:) = NaN;

ky2_cutoff(indf.i2,:,:) = NaN;
ky3_cutoff(indf.i3,:,:) = NaN;
ky4_cutoff(indf.i4,:,:) = NaN;

kym2(indf.i2,:) = NaN;
kym3(indf.i3,:) = NaN;
kym4(indf.i4,:) = NaN;

Sky2_wb(indf.i2,:,:) = NaN;
Sky3_wb(indf.i3,:,:) = NaN;
Sky4_wb(indf.i4,:,:) = NaN;

mag_wb2(indf.i2,:) = NaN;
mag_wb3(indf.i3,:) = NaN;
mag_wb4(indf.i4,:) = NaN;

wkym2(indf.i2,:) = NaN;
wkym3(indf.i3,:) = NaN;
wkym4(indf.i4,:) = NaN;

ky2_wb(indf.i2,:,:) = NaN;
ky3_wb(indf.i3,:,:) = NaN;
ky4_wb(indf.i4,:,:) = NaN;

wfit2.So(indf.i2,:) = NaN;wfit2.intw(indf.i2,:)=NaN;
wfit3.So(indf.i3,:) = NaN;wfit3.intw(indf.i3,:)=NaN;
wfit4.So(indf.i4,:) = NaN;wfit4.intw(indf.i4,:)=NaN;



head = 'get_weibul_take_out_bad.m';
save('/data1/bliu/data/edited_Wall','A','SS','head')

save('/data1/bliu/data/Sky_withWBfit_qced','ky_full','Sky2_cutoff_runmean','Sky3_cutoff_runmean','Sky4_cutoff_runmean','magspec2','magspec3','magspec4','Skym2','Skym3','Skym4','ky2_cutoff','ky3_cutoff',"ky4_cutoff",'kym2','kym3','kym4', ...
    'Sky2_wb','Sky3_wb','Sky4_wb','mag_wb2','mag_wb3','mag_wb4','wkym2','wkym3','wkym4','ky2_wb','ky3_wb','ky4_wb','wfit2','wfit3','wfit4','head')