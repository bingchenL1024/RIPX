%This code estimates the curl of wave breaking force from the 1.5 times
%the surf_zone width to the shoreline for 120 simulations
clear all
addpath(genpath('/data1/nkumar/RIPX/M_Files/funwavec_suite'));

%Load Ata's RIPX bathymtery guide
load RIPX_bath_guide

%Create a loop for number of runs
NRUNS = 120;
NT    = 3000:7800;
L     = length(NT);
DX    = 1;
DY    = 1.00;
LY    = 1200;

for i=1:1:NRUNS
    if i<10
        fname = sprintf('/data1/ssuanda/RIPX/run_%d/snap_%04d_',i,i);
        dname = sprintf('/data1/ssuanda/RIPX/run_%d/',i);
    else
        fname = sprintf('/data/ssuanda2/RIPX/run_%d/snap_%04d_',i,i);
        dname = sprintf('/data/ssuanda2/RIPX/run_%d/',i);
    end
    
    %Get water depth
    str0 = sprintf('%d',i);
    eval(['fname0=','runlist.bath.bath_',str0,';']);
    str0   = [dname,fname0];
    eval(['Depth=load(','''',str0,'''',');'])
    Depth  = repmat(Depth',1,1200);
    
    %Get cross-shore locations where we want to save stuff
    LX  = twocol(i,1):twocol(i,2);
    
    %Loop from 3001 to 7800 seconds
    for j=1:1:L
        str = sprintf('%04d',NT(j));
        %Get ETA
        eval(['fname1=','''',fname,'eta_snap_l2_',str,'.mat','''',';']);
        eval(['eta = load(','''',fname1,'''',');']);
        %Get EDDY
        eval(['fname1=','''',fname,'eddy_snap_l2_',str,'.mat','''',';']);
        eval(['eddy = load(','''',fname1,'''',');']);
        %Get U
        eval(['fname1=','''',fname,'u_snap_l2_',str,'.mat','''',';']);
        eval(['u = load(','''',fname1,'''',');']);
        %Get V
        eval(['fname1=','''',fname,'v_snap_l2_',str,'.mat','''',';']);
        eval(['v = load(','''',fname1,'''',');']);
        %Get VORT
        eval(['fname1=','''',fname,'vort_snap_l2_',str,'.mat','''',';']);
        eval(['vort = load(','''',fname1,'''',');']);
        
        [D_X, D_Y]       = funwaveC_Fbreaking(u.A,v.A,Depth+eta.A,eddy.A,DX,DY);
        DUM1             = funwaveC_curl(D_X, D_Y, DX, DY);
        [DUM2,~,~,~,~,~] = get_vel_decomposition(D_X',D_Y',DX,DY);
        psi(:,:,j)       = DUM2(:,LX)';
        curlF(:,:,j)     = DUM1(LX,:);
        disp(sprintf('Doing %d of %d',j,L));
        clear eta eddy u v vort D_X D_Y DUM1 DUM2 
    end
    disp(sprintf('Doing %d of %d',i,NRUNS));
    str0  = sprintf('%04d',i);
    sname = ['RIPX_psi_curlF_',str0,'.mat'];
    eval(['save -v7.3 ',sname,' psi',' curlF'])
    clear psi curlF
end
