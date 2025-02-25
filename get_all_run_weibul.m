% Bingchen Liu June 30, 2024
% modified from Ata's original code 


% Go through and save all params across the surfzone for all runs:
warning off
addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))

%load /data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat % now S(1--120) struct is in memory
load ('/data1/bliu/data/SS_raw.mat')
%%
% This script is basically run 3 times for each of the runs collected by
% slope. Loops through every 5m through the nearshore (includes locations outside of the surfzone)

slp = 2;
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, slp*.01);
for N = 1:24; % original one
%for N = 17:24; % BL modified for testing, change it back to avoid bug
ip4 = p4(N);
SS4 = S(ip4);  % all model variables are now in SS4

% This loads the spectra of breaking wave forcing
P4 = load_RIPX_curlF_Sky( ip4 );
%Grab the spectra:
Sky4 = P4.Sckm;
kyy = P4.fck/(2*pi);
xlim2 = [0.001 0.2];
% Need to subset to every 5m in surf:
[C,IA,IB] = intersect(SS4.X,SS4.X2); %29 is the cross-shore dimension 
Sky4 = Sky4(IA,:);
[iss4] = get_3locs(SS4);

for xx = 1:length(IA);
[S4] = interpSky(Sky4(xx,:),kyy,91);
iS4s(xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));
% Do the fit:
[FP4] = rayleigh_fitAS2(S4.lky,S4.sSkyl,1.2);
lskil4(xx) = falk_skill(log10(FP4.line'),log10(S4.sSkyl));
wskil4(xx) = wilmot_skill(log10(FP4.line'),log10(S4.sSkyl));
df = diff(S4.lky);
intw(xx) = nansum(FP4.line.*[0 df]');
a = find(S4.sSkyl==max(S4.sSkyl));kym4(xx) = S4.lky(a(1));
a = find(FP4.line==max(FP4.line));wkym4(xx) = S4.lky(a(1));

% Scales for the fits:
% offshore cyclic wavenum:
om4 = v4(N,2).*(2.*pi);
h4(xx) = SS4.h2(IB(xx))';
Hs4(xx) = SS4.Hs(IB(xx))';
kl4(xx) = get_wavenum(om4,h4(xx));
Dw4(xx) = real(SS4.dECG(IB(xx)));
hb(xx) = SS4.hb;
Ir(xx) = SS4.Irbo;
beta = 0.02;
%disp('done')
end
kl4pm = (kl4./(2*pi));

X = SS4.X(IA);

A.beta = beta;
A.Ir = Ir;
A.inp.tp = (2*pi)./om4;
A.inp.slp = slp*.01;
A.inp.dsp = v4(N,4);
A.X3locs = SS4.X(iss4);
A.h3locs = SS4.h(iss4);
A.X = X';
A.h = h4;
A.Hs = Hs4;
A.hb = hb;
A.kp = kl4pm;
A.Dw = Dw4;
A.raw.is = iS4s;
A.raw.kyo = kym4;
A.web.is = intw;
A.web.kyo = wkym4;
A.web.lsk = lskil4;
A.web.wsk = wskil4;
clear wskil4 lskil4 wkym4 intw kym4 iS4s Dw4 kl4 kl4pm Hs4 h4 X iss4 

eval(sprintf('Al_s%1.0f.a%1.0f = A;',slp,N))
disp('done')
clearvars -except Al_s2 S slp p4 v4 N
end
%% Slope 0.03
slp = 3;
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, slp*.01);

for N = 1:24; 
ip4 = p4(N);
SS4 = S(ip4);  % all model variables are now in SS4

% This loads the spectra of breaking wave forcing
P4 = load_RIPX_curlF_Sky( ip4 );
%Grab the spectra:
Sky4 = P4.Sckm;
kyy = P4.fck/(2*pi);
xlim2 = [0.001 0.2];
% Need to subset to every 5m in surf:
[C,IA,IB] = intersect(SS4.X,SS4.X2);
Sky4 = Sky4(IA,:);
[iss4] = get_3locs(SS4);

for xx = 1:length(IA);
[S4] = interpSky(Sky4(xx,:),kyy,91);
iS4s(xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));
% Do the fit:
[FP4] = rayleigh_fitAS2(S4.lky,S4.sSkyl,1.2);
lskil4(xx) = falk_skill(log10(FP4.line'),log10(S4.sSkyl));
wskil4(xx) = wilmot_skill(log10(FP4.line'),log10(S4.sSkyl));
df = diff(S4.lky);
intw(xx) = nansum(FP4.line.*[0 df]');
a = find(S4.sSkyl==max(S4.sSkyl));kym4(xx) = S4.lky(a(1));
a = find(FP4.line==max(FP4.line));wkym4(xx) = S4.lky(a(1));

% Scales for the fits:
% offshore cyclic wavenum:
om4 = v4(N,2).*(2.*pi);
h4(xx) = SS4.h2(IB(xx))';
Hs4(xx) = SS4.Hs(IB(xx))';
kl4(xx) = get_wavenum(om4,h4(xx));
Dw4(xx) = real(SS4.dECG(IB(xx)));
hb(xx) = SS4.hb;
Ir(xx) = SS4.Irbo;
beta = 0.03;
%disp('done')
end
kl4pm = (kl4./(2*pi));

X = SS4.X(IA);

A.beta = beta;
A.Ir = Ir;
A.inp.tp = (2*pi)./om4;
A.inp.slp = slp*.01;
A.inp.dsp = v4(N,4);
A.X3locs = SS4.X(iss4);
A.h3locs = SS4.h(iss4);
A.X = X';
A.h = h4;
A.Hs = Hs4;
A.hb= hb;
A.kp = kl4pm;
A.Dw = Dw4;
A.raw.is = iS4s;
A.raw.kyo = kym4;
A.web.is = intw;
A.web.kyo = wkym4;
A.web.lsk = lskil4;
A.web.wsk = wskil4;
clear wskil4 lskil4 wkym4 intw kym4 iS4s Dw4 kl4 kl4pm Hs4 h4 X iss4 

eval(sprintf('Al_s%1.0f.a%1.0f = A;',slp,N))
disp('done')
clearvars -except Al_s2 Al_s3 S slp p4 v4 N
end
%%
% Slope 0.04
slp = 4;
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, slp*.01);

for N = 1:24; 
ip4 = p4(N);
SS4 = S(ip4);  % all model variables are now in SS4

% This loads the spectra of breaking wave forcing
P4 = load_RIPX_curlF_Sky( ip4 );
%Grab the spectra:
Sky4 = P4.Sckm;
kyy = P4.fck/(2*pi);
xlim2 = [0.001 0.2];
% Need to subset to every 5m in surf:
[C,IA,IB] = intersect(SS4.X,SS4.X2);
Sky4 = Sky4(IA,:);
[iss4] = get_3locs(SS4);

for xx = 1:length(IA);
[S4] = interpSky(Sky4(xx,:),kyy,91);
iS4s(xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));
% Do the fit:
[FP4] = rayleigh_fitAS2(S4.lky,S4.sSkyl,1.2);
lskil4(xx) = falk_skill(log10(FP4.line'),log10(S4.sSkyl));
wskil4(xx) = wilmot_skill(log10(FP4.line'),log10(S4.sSkyl));
df = diff(S4.lky);
intw(xx) = nansum(FP4.line.*[0 df]');
a = find(S4.sSkyl==max(S4.sSkyl));kym4(xx) = S4.lky(a(1));
a = find(FP4.line==max(FP4.line));wkym4(xx) = S4.lky(a(1));

% Scales for the fits:
% offshore cyclic wavenum:
om4 = v4(N,2).*(2.*pi);
h4(xx) = SS4.h2(IB(xx))';
Hs4(xx) = SS4.Hs(IB(xx))';
kl4(xx) = get_wavenum(om4,h4(xx));
Dw4(xx) = real(SS4.dECG(IB(xx)));
hb(xx) = SS4.hb;
Ir(xx) = SS4.Irbo;
beta = 0.04;
%disp('done')
end
kl4pm = (kl4./(2*pi));

X = SS4.X(IA);

A.beta = beta;
A.Ir = Ir;
A.inp.tp = (2*pi)./om4;
A.inp.slp = slp*.01;
A.inp.dsp = v4(N,4);
A.X3locs = SS4.X(iss4);
A.h3locs = SS4.h(iss4);
A.X = X';
A.h = h4;
A.Hs = Hs4;
A.hb = hb;
A.kp = kl4pm;
A.Dw = Dw4;
A.raw.is = iS4s;
A.raw.kyo = kym4;
A.web.is = intw;
A.web.kyo = wkym4;
A.web.lsk = lskil4;
A.web.wsk = wskil4;
clear wskil4 lskil4 wkym4 intw kym4 iS4s Dw4 kl4 kl4pm Hs4 h4 X iss4 

eval(sprintf('Al_s%1.0f.a%1.0f = A;',slp,N))
disp('done')
clearvars -except Al_s2 Al_s3 Al_s4 S slp p4 v4 N
end

head = 'All profiles, different fit ranges. see get_all_run_weibul.m';
save('/data1/bliu/data/Profs_fit_to0p2_BL','Al_s4','Al_s3','Al_s2','head')
