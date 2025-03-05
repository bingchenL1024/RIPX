% Bingchen Liu Nov 21, 2024
% This code get the nondimensional spectra (full)
% It normalize the spectra by its variance and ky by peak kym


clearvars -except expo


%load('/data1/bliu/data/Sky_withWBfit_qced') %original default expo = 1.2



%expo=1.25;
expo_name_temp = num2str(expo);
expo_name  = [expo_name_temp(1),'p',expo_name_temp(3:end)];

load(['/data1/bliu/data/Sky_withWBfit_qced_',expo_name])

for N = 1:24
    for xx=1:3
        % Sky2_nond(N,xx,:) = Sky2_cutoff_runmean(N,xx,:)./magspec2(N,xx);
        % Sky3_nond(N,xx,:) = Sky3_cutoff_runmean(N,xx,:)./magspec3(N,xx);
        % Sky4_nond(N,xx,:) = Sky4_cutoff_runmean(N,xx,:)./magspec4(N,xx);
        
        % Sky2_nond_maxsky(N,xx,:) = Sky2_cutoff_runmean(N,xx,:)./Skym2(N,xx);
        % Sky3_nond_maxsky(N,xx,:) = Sky3_cutoff_runmean(N,xx,:)./Skym3(N,xx);
        % Sky4_nond_maxsky(N,xx,:) = Sky4_cutoff_runmean(N,xx,:)./Skym4(N,xx);
        % 
        % Sky2_nond_kym(N,xx,:) = Sky2_cutoff_runmean(N,xx,:).*kym2(N,xx)./magspec2(N,xx);
        % Sky3_nond_kym(N,xx,:) = Sky3_cutoff_runmean(N,xx,:).*kym3(N,xx)./magspec3(N,xx);
        % Sky4_nond_kym(N,xx,:) = Sky4_cutoff_runmean(N,xx,:).*kym4(N,xx)./magspec4(N,xx);
        % 
        % Sky2_nond_kym_S0(N,xx,:) = Sky2_cutoff_runmean(N,xx,:).*kym2(N,xx)./wfit2.So(N,xx);
        % Sky3_nond_kym_S0(N,xx,:) = Sky2_cutoff_runmean(N,xx,:).*kym3(N,xx)./wfit3.So(N,xx); 
        % Sky4_nond_kym_S0(N,xx,:) = Sky2_cutoff_runmean(N,xx,:).*kym4(N,xx)./wfit4.So(N,xx);

        Sky2_nond_kym_S0(N,xx,:) = Sky2_cutoff_runmean(N,xx,:).*wkym2(N,xx)./wfit2.So(N,xx);
        Sky3_nond_kym_S0(N,xx,:) = Sky3_cutoff_runmean(N,xx,:).*wkym3(N,xx)./wfit3.So(N,xx); 
        Sky4_nond_kym_S0(N,xx,:) = Sky4_cutoff_runmean(N,xx,:).*wkym4(N,xx)./wfit4.So(N,xx);

        ky2_nond(N,xx,:) = ky2_cutoff(N,xx,:)/kym2(N,xx);
        ky3_nond(N,xx,:) = ky3_cutoff(N,xx,:)/kym3(N,xx);
        ky4_nond(N,xx,:) = ky4_cutoff(N,xx,:)/kym4(N,xx);

        % Sky2_nond_kym_wb(N,xx,:) = Sky2_wb(N,xx,:).*wkym2(N,xx)./mag_wb2(N,xx);
        % Sky3_nond_kym_wb(N,xx,:) = Sky3_wb(N,xx,:).*wkym3(N,xx)./mag_wb3(N,xx);
        % Sky4_nond_kym_wb(N,xx,:) = Sky4_wb(N,xx,:).*wkym4(N,xx)./mag_wb4(N,xx);
        % 
        % Sky2_nond_kym_wb_S0(N,xx,:) = Sky2_wb(N,xx,:).*wkym2(N,xx)./wfit2.So(N,xx);
        % Sky3_nond_kym_wb_S0(N,xx,:) = Sky3_wb(N,xx,:).*wkym3(N,xx)./wfit3.So(N,xx);
        % Sky4_nond_kym_wb_S0(N,xx,:) = Sky4_wb(N,xx,:).*wkym4(N,xx)./wfit4.So(N,xx);



        % ky2_nond_wb(N,xx,:) = ky2_wb(N,xx,:)/wkym2(N,xx);
        % ky3_nond_wb(N,xx,:) = ky3_wb(N,xx,:)/wkym3(N,xx);
        % ky4_nond_wb(N,xx,:) = ky4_wb(N,xx,:)/wkym4(N,xx);

    end 
end 

save(['/data1/bliu/data/Sky_nond_',expo_name])


%% testing and debugging 
% for xloc=1:3
% test_spec = squeeze(Sky2_nond(N,:,:));
% test_freq = squeeze(ky2_cutoff(N,:,:));
% nansum(test_spec(xloc,:)).*diff(test_freq(xloc,1:2))
% end 
