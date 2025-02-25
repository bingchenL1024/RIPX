% Calculate the scaled So at the three locations:
warning off
addpath('/home/ffeddersen/matlab')
addpath(genpath('/data1/nkumar/RIPX/M_Files'))

load /data1/nkumar/RIPX/M_Files/Important_Files/qa_qc_RIPX_NK.mat % now S(1--120) struct is in memory

% this returns the RIPX run # for the requested parameters.
[p2,v2] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.02);
[p3,v3] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.03);
[p4,v4] = get_sampled_ripx_spread_index(0.1, 2, 0.01, 0.5, 0.04);

% Just grab the offshore cyclic omega:
om2 = v2(:,2).*(2.*pi);
om3 = v3(:,2).*(2.*pi);
om4 = v4(:,2).*(2.*pi);

for N=1:length(p4)
    ip2 = p2(N);
    ip3 = p3(N);
    ip4 = p4(N);
SS2 = S(ip2);  % all model variables are now in SS2
SS3 = S(ip3);  % all model variables are now in SS3
SS4 = S(ip4);  % all model variables are now in SS4

% Get the indexes for a 3 different locations within surfzone 
[iss2] = get_3locs(SS2,2);
[iss3] = get_3locs(SS3,2);
[iss4] = get_3locs(SS4,2);

%Grab the water depth:
h2(N,:) = SS2.h2(iss2)';
h3(N,:) = SS3.h2(iss3)';
h4(N,:) = SS4.h2(iss4)';
%Grab the wavenumber:
kl2(N,:) = get_wavenum(om2(N),h2(N,:));
kl3(N,:) = get_wavenum(om3(N),h3(N,:));
kl4(N,:) = get_wavenum(om4(N),h4(N,:));
% Grab the dissipation:
Dw2(N,:) = SS2.dECG(iss2);
Dw3(N,:) = SS3.dECG(iss3);
Dw4(N,:) = SS4.dECG(iss4);
end

% So, that is in rad/m: this one is in 1/m:
kl2pm = (kl2./(2*pi));
kl3pm = (kl3./(2*pi));
kl4pm = (kl4./(2*pi));

save('Syscale','kl2','kl3','kl4','kl2pm','kl3pm','kl4pm','h2','h3','h4','Dw2','Dw3','Dw4')