% This script does the Sky truncating, interpolating, and Weibul fitting for 3 locations in surfzone
% A .mat file is saved at the end with integrated spectra, fits, and skill scores
%
% Step 1: Get all the Weibul fits
warning off

addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))

load /data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat % now S(1--120) struct is in memory

% this returns the RIPX run # for the requested parameters.
[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

for N=1:length(p4)
    ip2 = p2(N); %index number with given bath 
    ip3 = p3(N);
    ip4 = p4(N);
SS2 = S(ip2);  % all model variables are now in SS2
SS3 = S(ip3);  % all model variables are now in SS3
SS4 = S(ip4);  % all model variables are now in SS4

% This loads the spectra of breaking wave forcing
P2 = load_RIPX_curlF_Sky( ip2 );
P3 = load_RIPX_curlF_Sky( ip3 );
P4 = load_RIPX_curlF_Sky( ip4 );

% Get the indexes for a 3 different locations within surfzone 
[iss2] = get_3locs(SS2);
[iss3] = get_3locs(SS3);
[iss4] = get_3locs(SS4);

%Grab the spectra:
Sky2 = P2.Sckm(iss2,:);Sky3 = P3.Sckm(iss3,:);Sky4 = P4.Sckm(iss4,:);
kyy = P2.fck/(2*pi);
xlim2 = [0.0008 0.2];

% For the three locations:
    for xx = 1:3;
        % 1. Truncation and logarithmic interpolation:
[S2] = interpSky(Sky2(xx,:),kyy,91); % 91 is width of log-bins
[S3] = interpSky(Sky3(xx,:),kyy,91);
[S4] = interpSky(Sky4(xx,:),kyy,91);

% 2. Get the total variance in Sky and its integration out to 0.2:
iS2tot(N,xx) = (sum(Sky2(xx,:)).*diff(kyy(1:2))); % integrate the whole thing
iS3tot(N,xx) = (sum(Sky3(xx,:)).*diff(kyy(1:2)));
iS4tot(N,xx) = (sum(Sky4(xx,:)).*diff(kyy(1:2)));
iS2s(N,xx) = (nansum(S2.Skyi).*diff(S2.kyy(1:2))); % only to 0.2
iS3s(N,xx) = (nansum(S3.Skyi).*diff(S3.kyy(1:2)));
iS4s(N,xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));

% Do the Weibul fits:
expo = 1.2; % Choose the exponent
[FP2] = rayleigh_fitAS2(S2.lky,S2.sSkyl,expo);
[FP3] = rayleigh_fitAS2(S3.lky,S3.sSkyl,expo);
[FP4] = rayleigh_fitAS2(S4.lky,S4.sSkyl,expo);

% 2 metrics for the skill of the fit:
lskil2(N,xx) = falk_skill(log10(FP2.line'),log10(S2.sSkyl));
lskil3(N,xx) = falk_skill(log10(FP3.line'),log10(S3.sSkyl));
lskil4(N,xx) = falk_skill(log10(FP4.line'),log10(S4.sSkyl));
wskil2(N,xx) = wilmot_skill(log10(FP2.line'),log10(S2.sSkyl));
wskil3(N,xx) = wilmot_skill(log10(FP3.line'),log10(S3.sSkyl));
wskil4(N,xx) = wilmot_skill(log10(FP4.line'),log10(S4.sSkyl));

% Grab the maximum value and kyy value:
a = find(S2.sSkyl==max(S2.sSkyl)); Skym2(N,xx) = S2.sSkyl(a); kym2(N,xx) = S2.lky(a);
a = find(S3.sSkyl==max(S3.sSkyl)); Skym3(N,xx) = S3.sSkyl(a); kym3(N,xx) = S3.lky(a);
a = find(S4.sSkyl==max(S4.sSkyl)); Skym4(N,xx) = S4.sSkyl(a); kym4(N,xx) = S4.lky(a);
a = find(FP2.line==max(FP2.line)); wSkym2(N,xx) = FP2.line(a); wkym2(N,xx) = S2.lky(a);
a = find(FP3.line==max(FP3.line)); wSkym3(N,xx) = FP3.line(a); wkym3(N,xx) = S3.lky(a);
a = find(FP4.line==max(FP4.line)); wSkym4(N,xx) = FP4.line(a); wkym4(N,xx) = S4.lky(a);

% Also save the individual fit params a and b:
wfit2.a(N,xx) = FP2.a; wfit2.b(N,xx) = FP2.b;
wfit3.a(N,xx) = FP3.a; wfit3.b(N,xx) = FP3.b;
wfit4.a(N,xx) = FP4.a; wfit4.b(N,xx) = FP4.b;

% Integrate the Weibul:
df = diff(S2.lky);
wfit2.intw(N,xx) = nansum(FP2.line.*[0 df]');
df = diff(S3.lky);
wfit3.intw(N,xx) = nansum(FP3.line.*[0 df]');
df = diff(S4.lky);
wfit4.intw(N,xx) = nansum(FP4.line.*[0 df]');
% If the whole Weibul part integrates to 1, then the variance should be:
wfit2.So(N,xx) = FP2.a./expo;
wfit3.So(N,xx) = FP3.a./expo;
wfit4.So(N,xx) = FP4.a./expo;

    end
    disp('done')
end
save('WF_1p2','lskil2','lskil3','lskil4','wskil2','wskil3','wskil4','Skym2','Skym3','Skym4','kym2','kym3','kym4','wSkym2','wSkym3','wSkym4','wkym2','wkym3','wkym4','wfit2','wfit3','wfit4','iS2tot','iS3tot','iS4tot','iS2s','iS3s','iS4s')