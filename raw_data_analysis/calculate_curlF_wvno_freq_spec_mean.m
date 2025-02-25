clear all
addpath(genpath('/home/n2kumar/MAT_Mytoolbox/jlab/'));
NT      = 120;
NS      = 1;
dirname = '/data1/nkumar/RIPX/M_Files/RIPX_curlF_wvno_spec/';
dirsname= '/data1/nkumar/RIPX/M_Files/RIPX_mean_curlF_kf_mspec/';
dtime   = 3; %average every 3 time stamps

for i=NS:1:NT
    str0         = sprintf('%04d',i);
    fname1       = [dirname,'RIPX_curlF_wvno_spec_',str0,'.mat'];
    eval(['Fname1=matfile(','''',fname1,'''',');']);
    
    Sname        = [dirsname,'RIPX_mean_curlF_kf_mspec_',str0,'.mat'];
    
    disp('Getting mean of wave no spectra');
    Sckm         = nanmean(Fname1.Sck(:,1:dtime:end,:),2);% take the mean in time for the spectra
    Sckm         = squeeze(Sckm(:,1,:));
    fck          = Fname1.fck;
    eval(['save ',Sname,' Sckm',' fck']);
    disp(i)
end
