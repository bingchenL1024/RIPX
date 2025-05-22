% Bingchen Liu June 30, 2024
% modified to include 10 different surfzone locations  

% Get SS from  basic data from S.mat + incoming wave spectra (include 3 locations)
clear
warning off
addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))

load /data1/bliu/data/SS_raw.mat % now S(1--120) struct is in memory
load('/data1/bliu/data/RIPX_bath_guide.mat')
% this returns the RIPX run # for the requested parameters.
[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02); %p2 is runnumber; v2=[Hs3 fp3 m3 spread3];
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

% Just grab the offshore cyclic omega:
om2 = v2(:,2).*(2.*pi);
om3 = v3(:,2).*(2.*pi);
om4 = v4(:,2).*(2.*pi);

for N=1:length(p4) %N = 1:24(different run with the same bath)
    ip2 = p2(N); %ip2 is index of runnumber of the Nth in the same bath;
    ip3 = p3(N);
    ip4 = p4(N);
SS2 = S(ip2);  % all model variables are now in SS2, get a slice of S with diff bath
SS3 = S(ip3);  % all model variables are now in SS3
SS4 = S(ip4);  % all model variables are now in SS4


bath2= AA.A(ip2,1);
bath3= AA.A(ip3,1);
bath4= AA.A(ip4,1);

% Get the indexes for a 3 different locations within surfzone 
[iss2] = get_10locs(SS2,2); 
[iss3] = get_10locs(SS3,2);
[iss4] = get_10locs(SS4,2);

[iss2_1] = get_10locs(SS2,1); 
[iss3_1] = get_10locs(SS3,1);
[iss4_1] = get_10locs(SS4,1);
%Grab the water depth:
SS.h2(N,:) = SS2.h2(iss2)';
SS.h3(N,:) = SS3.h2(iss3)';
SS.h4(N,:) = SS4.h2(iss4)';
% Grab water depth at breaking
SS.hb2(N,1) = SS2.hb;
SS.hb3(N,1) = SS3.hb;
SS.hb4(N,1) = SS4.hb;

%Grab the wavenumber:
SS.kl2(N,:) = get_wavenum(om2(N),SS.h2(N,:)); % use omega at breaking pt 
SS.kl3(N,:) = get_wavenum(om3(N),SS.h3(N,:));
SS.kl4(N,:) = get_wavenum(om4(N),SS.h4(N,:));
% Grab the dissipation:
SS.Dw2(N,:) = SS2.dECG(iss2);
SS.Dw3(N,:) = SS3.dECG(iss3);
SS.Dw4(N,:) = SS4.dECG(iss4);

% Grab the curlF_std
SS.curlFstd2(N,:) = SS2.curlF_std(iss2_1);
SS.curlFstd3(N,:) = SS3.curlF_std(iss3_1);
SS.curlFstd4(N,:) = SS4.curlF_std(iss4_1);

% grab wave steepness at breaking
SS.stpb2(N,:) = SS2.stpb;
SS.stpb3(N,:) = SS3.stpb;
SS.stpb4(N,:) = SS4.stpb;


%Grab Iribarren number
SS.Ir2(N,:) = SS2.Irbo;
SS.Ir3(N,:) = SS3.Irbo;
SS.Ir4(N,:) = SS4.Irbo;

SS.Irb2(N,:) = 0.02/(SS2.stpb)^0.5;
SS.Irb3(N,:) = 0.03/(SS3.stpb)^0.5;
SS.Irb4(N,:) = 0.04/(SS4.stpb)^0.5;


% Grab the directional spread:
ds.sigtb2(N) = SS2.sigma_b;
ds.sigtb3(N) = SS3.sigma_b;
ds.sigtb4(N) = SS4.sigma_b;

% Grab beach slope/bath (beta)
SS.beta2(N,1) = bath2;
SS.beta3(N,1) = bath3;
SS.beta4(N,1) = bath4;

end

% peak wave number of incoming wave spectra: this one is in 1/m:
SS.kl2pm = (SS.kl2./(2*pi));
SS.kl3pm = (SS.kl3./(2*pi));
SS.kl4pm = (SS.kl4./(2*pi));

README = "Data generated from 'get_SS_ds_10loc.m'";
README = [README "Include basic data from S.mat(include 10 locations) + directional spread (same for all 10 locs)"];
README = splitlines(README)';

save('/data1/bliu/data/SS_ds_10loc','README','SS','ds')