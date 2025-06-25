% Modified by Bingchen Liu Jan 23, 2025
% modified fomr get_Sky_weibul_fits_log.mat' --> that saves output 'WF_1p2.mat'
% current output is 'Sky_WBfit_expo_1p2'
% Include saving the raw spectra as well as the truncated spectra 
% Note the original code does not save any raw spectra info, it only saves
% fit parameter info
% used in 'master_analysis_bl.m'
% compare to 'get_Sky_weibul_fits_log_BL_testexp', this code include expo
% inside the script --> don't need outside input --> it can run indepdently

% This script does the Sky truncating, interpolating, and Weibul fitting for 3 locations in surfzone
% A .mat file is saved at the end with integrated spectra, fits, and skill scores
%
% Step 1: Get all the Weibul fits
warning off
%clearvars -except expo badfit_num diffexpo rmse_diffexp rmse_diffexp_log rmse_diffexp_weight

addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))

load /data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat % now S(1--120) struct is in memory



%expo = 1.25; % Choose the exponent



% this returns the RIPX run # for the requested parameters.
[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

for N=1:length(p4)
    ip2 = p2(N);
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
[iss2] = get_10locs(SS2);
[iss3] = get_10locs(SS3);
[iss4] = get_10locs(SS4);

%Grab the spectra:
Sky2 = P2.Sckm(iss2,:);Sky3 = P3.Sckm(iss3,:);Sky4 = P4.Sckm(iss4,:);
kyy = P2.fck/(2*pi); % kyy (CPM)
xlim2 = [0.0008 0.2];

loc_num = length(iss2);
% For the three locations:
    for xx = 1:loc_num
        % 1. Truncation and logarithmic interpolation:
[S2] = interpSky(Sky2(xx,:),kyy,91); % 91 is width of log-bins for running mean to smooth the spectra
[S3] = interpSky(Sky3(xx,:),kyy,91);
[S4] = interpSky(Sky4(xx,:),kyy,91);

% 2. Get the total variance in Sky and its integration out to 0.2:
iS2tot(N,xx) = (sum(Sky2(xx,3:end)).*diff(kyy(1:2))); % integrate the whole thing
iS3tot(N,xx) = (sum(Sky3(xx,3:end)).*diff(kyy(1:2)));
iS4tot(N,xx) = (sum(Sky4(xx,3:end)).*diff(kyy(1:2)));
iS2s(N,xx) = (nansum(S2.Skyi).*diff(S2.kyy(1:2))); % only to 0.2
iS3s(N,xx) = (nansum(S3.Skyi).*diff(S3.kyy(1:2)));
iS4s(N,xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));

magspec2_temp(xx) = (nansum(S2.Skyi).*diff(S2.kyy(1:2))); % just for saving in different format for plot purpose
magspec3_temp(xx) = (nansum(S3.Skyi).*diff(S3.kyy(1:2))); 
magspec4_temp(xx) = (nansum(S4.Skyi).*diff(S4.kyy(1:2)));

% magspec2_temp(xx) = (sum(Sky2(xx,:)).*diff(kyy(1:2))); % just for saving in different format for plot purpose
% magspec3_temp(xx) = (sum(Sky3(xx,:)).*diff(kyy(1:2))); 
% magspec4_temp(xx) = (sum(Sky4(xx,:)).*diff(kyy(1:2))); 


Sky2_cutoff_temp(xx,:) = S2.sSkyi;
Sky3_cutoff_temp(xx,:) = S3.sSkyi;
Sky4_cutoff_temp(xx,:) = S4.sSkyi;

kyy2_cutoff_temp(xx,:) = S2.kyy;
kyy3_cutoff_temp(xx,:) = S3.kyy;
kyy4_cutoff_temp(xx,:) = S4.kyy;

kyy2_wb_temp(xx,:) = S2.lky;
kyy3_wb_temp(xx,:) = S3.lky;
kyy4_wb_temp(xx,:) = S4.lky;


ky_full=kyy;
% ===================>>>>>>>>>>>> Do the Weibul fits<<<<<<<<==============
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

% rmse_error2(N,xx) = rmse(FP2.line',S2.sSkyl);
% rmse_error3(N,xx) = rmse(FP3.line',S3.sSkyl);
% rmse_error4(N,xx) = rmse(FP4.line',S4.sSkyl);
% 
% rmse_error_log2(N,xx) = rmse(log(FP2.line'),log(S2.sSkyl));
% rmse_error_log3(N,xx) = rmse(log(FP3.line'),log(S3.sSkyl));
% rmse_error_log4(N,xx) = rmse(log(FP4.line'),log(S4.sSkyl));
% 
% rmse_weight2(N,xx) = iS2s(N,xx);
% rmse_weight3(N,xx) = iS3s(N,xx);
% rmse_weight4(N,xx) = iS4s(N,xx);


% Grab the maximum value and kyy value:
a = find(S2.sSkyl==max(S2.sSkyl));Skym2(N,xx) = S2.sSkyl(a);kym2(N,xx) = S2.lky(a);
a = find(S3.sSkyl==max(S3.sSkyl));Skym3(N,xx) = S3.sSkyl(a);kym3(N,xx) = S3.lky(a);
a = find(S4.sSkyl==max(S4.sSkyl));Skym4(N,xx) = S4.sSkyl(a);kym4(N,xx) = S4.lky(a);
a = find(FP2.line==max(FP2.line));wSkym2(N,xx) = FP2.line(a);wkym2(N,xx) = S2.lky(a);
a = find(FP3.line==max(FP3.line));wSkym3(N,xx) = FP3.line(a);wkym3(N,xx) = S3.lky(a);
a = find(FP4.line==max(FP4.line));wSkym4(N,xx) = FP4.line(a);wkym4(N,xx) = S4.lky(a);

% Also save the individual fit params a and b:
wfit2.a(N,xx) = FP2.a;wfit2.b(N,xx) = FP2.b;
wfit3.a(N,xx) = FP3.a;wfit3.b(N,xx) = FP3.b;
wfit4.a(N,xx) = FP4.a;wfit4.b(N,xx) = FP4.b;

% Integrate the Weibul:
df = diff(S2.lky);
wfit2.intw(N,xx) = nansum(FP2.line(219:end).*[df(218:end)]');
df = diff(S3.lky);
wfit3.intw(N,xx) = nansum(FP3.line(219:end).*[df(218:end)]');
df = diff(S4.lky);
wfit4.intw(N,xx) = nansum(FP4.line(219:end).*[df(218:end)]');

mag_wb2_temp(xx) = nansum(FP2.line.*[0 df]');
mag_wb3_temp(xx) = nansum(FP3.line.*[0 df]');
mag_wb4_temp(xx) = nansum(FP4.line.*[0 df]');
% If the whole Weibul part integrates to 1, then the variance should be:
wfit2.So(N,xx) = FP2.a./expo;
wfit3.So(N,xx) = FP3.a./expo;
wfit4.So(N,xx) = FP4.a./expo;

Sky2_wb_temp(xx,:) = FP2.line;
Sky3_wb_temp(xx,:) = FP3.line;
Sky4_wb_temp(xx,:) = FP4.line;

    end % 3 locations

Sky2_raw(N,:,:) = Sky2;
Sky3_raw(N,:,:) = Sky3;
Sky4_raw(N,:,:) = Sky4;


Sky2_cutoff_runmean(N,:,:) = Sky2_cutoff_temp;
Sky3_cutoff_runmean(N,:,:) = Sky3_cutoff_temp;
Sky4_cutoff_runmean(N,:,:) = Sky4_cutoff_temp;

magspec2(N,:,:)= magspec2_temp';
magspec3(N,:,:)= magspec3_temp';
magspec4(N,:,:)= magspec4_temp';

ky2_cutoff(N,:,:) = kyy2_cutoff_temp;
ky3_cutoff(N,:,:) = kyy3_cutoff_temp;
ky4_cutoff(N,:,:) = kyy4_cutoff_temp;

Sky2_wb(N,:,:) = Sky2_wb_temp;
Sky3_wb(N,:,:) = Sky3_wb_temp;
Sky4_wb(N,:,:) = Sky4_wb_temp;

mag_wb2(N,:) = mag_wb2_temp;
mag_wb3(N,:) = mag_wb3_temp;
mag_wb4(N,:) = mag_wb4_temp;

ky2_wb(N,:,:) = kyy2_wb_temp;
ky3_wb(N,:,:) = kyy3_wb_temp;
ky4_wb(N,:,:) = kyy4_wb_temp;
    %disp('done')
end

% rmse_tot_temp = [rmse_error2(:);rmse_error3(:);rmse_error4(:)];
% rmse_weight_tot_temp = [rmse_weight2(:);rmse_weight3(:);rmse_weight4(:)];
% 
% global rmse_unweight rmse_log_unweight rmse_weight rmse_log_weight
% rmse_unweight = mean([rmse_error2(:);rmse_error3(:);rmse_error4(:)]);
% rmse_log_unweight = mean([rmse_error_log2(:);rmse_error_log3(:);rmse_error_log4(:)]);
% rmse_weight = mean(rmse_tot_temp.*rmse_weight_tot_temp,'all')/sum(rmse_weight_tot_temp,'all');
% 


%%
% expo_name_temp = num2str(expo);
% expo_name  = [expo_name_temp(1),'p',expo_name_temp(3:end)];
% include important fit parameter
save('/data1/bliu/data/Sky_WBfit_10loc','lskil2','lskil3','lskil4','wskil2','wskil3','wskil4','Skym2','Skym3','Skym4','kym2','kym3','kym4','wSkym2','wSkym3','wSkym4','wkym2','wkym3','wkym4','wfit2','wfit3','wfit4','iS2tot','iS3tot','iS4tot','iS2s','iS3s','iS4s','Sky2_raw','Sky3_raw','Sky4_raw','ky_full','Sky2_cutoff_runmean','Sky3_cutoff_runmean','Sky4_cutoff_runmean','magspec2','magspec3','magspec4','Skym2','Skym3','Skym4','ky2_cutoff','ky3_cutoff',"ky4_cutoff",'kym2','kym3','kym4', ...
    'Sky2_wb','Sky3_wb','Sky4_wb','mag_wb2','mag_wb3','mag_wb4','wkym2','wkym3','wkym4','ky2_wb','ky3_wb','ky4_wb','wfit2','wfit3','wfit4')

% include more parameters and exclude some parameter --> only for nond Sky
% plot 
% save(['/data1/bliu/data/Sky_withWBfit_expo_',expo_name],'Sky2_raw','Sky3_raw','Sky4_raw','ky_full','Sky2_cutoff_runmean','Sky3_cutoff_runmean','Sky4_cutoff_runmean','magspec2','magspec3','magspec4','Skym2','Skym3','Skym4','ky2_cutoff','ky3_cutoff',"ky4_cutoff",'kym2','kym3','kym4', ...
%     'Sky2_wb','Sky3_wb','Sky4_wb','mag_wb2','mag_wb3','mag_wb4','wkym2','wkym3','wkym4','ky2_wb','ky3_wb','ky4_wb','wfit2','wfit3','wfit4')