% Bingchen Liu Nov 19, 2024
% This code post process the raw data and low pass it 

clear all
addpath(genpath('/home/ffeddersen/matlab'))
runind=2;
tic
clearvars -except runind
fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',runind);
load(fname,'curlF')
%%
dim = size(curlF);
f_data = 1; %HZ
f_cutoff = 0.3; %HZ
t=-20:1:20;
f_lp = (1/(sqrt(2*pi)*3)).*exp(-0.5*(t./3).^2);

for i = 1:dim(1)
    for j = 1:dim(2)
        if i ==110&&j==100
            keyboard
        end 
        data_raw= squeeze(curlF(i,j,:));
        %data_lp = lowpass(data_raw,f_cutoff,f_data);
        data_lp = conv(data_raw,f_lp);
        curlF_lp(i,j,:) = data_lp; 
    end 
end 

save('/data1/bliu/data/raw/curlF_lowpass/run_2','curlF_lp','-v7.3')