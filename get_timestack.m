% Bingchen Liu Feb 26, 2025
% This code get the time stack for paper fig ploting 
clear
runind=14;

head = "data generated from get_timestack.m";
clearvars -except runind
fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_vort/RIPX_vort_%04d',runind);

load(fname,'vort')
%vort = permute(vort,[2 3 1]);
dim_vort= size(vort);


yind = 10;
tlim = 2000:2150;


vort_tstack = squeeze(vort(tlim,:,yind))';
save('/data1/bliu/data/tstack_vort_run14', "vort_tstack")


clearvars -except runind yind tlim
fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',runind);
load(fname,'curlF')
dim_curlFbr= size(curlF);

curlF_tstack = squeeze(curlF(:,yind,tlim));
save('/data1/bliu/data/tstack_curlF_run14', "curlF_tstack")