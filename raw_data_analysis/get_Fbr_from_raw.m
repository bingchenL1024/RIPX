% Bingchen Liu June 14, 2024
% This code obtain the alongshore averaged Fbr magnitude (include both x and y component)
% The input data is time averaged Fbr_x and Fbr_y with dim: cross-shore X
% along-shore

clear all
NT        = 120;
NS        = 1;
dirname   = '/data1/nkumar/RIPX/M_Files/RIPX_std_Fbr/';
dirsname  = '/data1/bliu/data/raw/Fbr_ymean_tstd/';
 
for i=NS:1:NT
    str0         = sprintf('%04d',i);
    fname1       = [dirname,'RIPX_std_Fbr_',str0,'.mat'];
    eval(['Fname1=matfile(','''',fname1,'''',');']);
    Sname        = [dirsname,'RIPX_Fbr_ymean_tstd',str0,'.mat'];
    Fbrx_ymean_tstd = squeeze(mean(Fname1.Fbrx_std,2,'omitnan'));
    Fbry_ymean_tstd = squeeze(mean(Fname1.Fbry_std,2,'omitnan'));
    Fbr_mag_ymean_tstd = vecnorm([Fbrx_ymean_tstd,Fbry_ymean_tstd],2,2);
    disp(i)
    eval(['save ',Sname,' Fbr_mag_ymean_tstd',' Fbrx_ymean_tstd',' Fbry_ymean_tstd'])
    clear Fbr_mag_ymean_tstd  Fbrx_ymean_tstd Fbry_ymean_tstd 
end