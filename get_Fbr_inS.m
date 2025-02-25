% Bingchen Liu June 14, 2024
% This code add Fbr info (mag, x component, y component) into general
% data info S in 'qa_qc_RIPX_NK' (or S (var name))

clear all 
load('/data1/bliu/data/qa_qc_RIPX_NK.mat')
NT        = 120;
NS        = 1;

for i=NS:1:NT
    dirname        = '/data1/bliu/data/raw/Fbr_ymean_tstd/'; %BL modified 
    str        = sprintf('%04d',i);
    Fname      = [dirname,'RIPX_Fbr_ymean_tstd',str,'.mat'];
    eval(['load ','''',Fname,'''']);
    
    S(i).Fbr_mag = Fbr_mag_ymean_tstd;
    S(i).Fbr_x = Fbrx_ymean_tstd;  
    S(i).Fbr_y = Fbry_ymean_tstd;
    
    clear Fbr_mag_ymean_tstd  Fbrx_ymean_tstd Fbry_ymean_tstd 

end 

README = "Data generated from 'get_Fbr_inS.m'";
README = [README "basically add Fbr to S (in 'qa_qc_RIPX_NK.mat')"];
README = splitlines(README)';
save('/data1/bliu/data/SS_raw.mat','S', 'README')
