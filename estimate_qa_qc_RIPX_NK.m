%This function estimates the energy flux, energy flux gradient, and many
%other terms for each simulation. 

%This code is modified by Bingchen Liu. The part that has been modified is
%labeled as BL edit 
addpath(genpath('/home/n2kumar/MTOOLS/'));
addpath(genpath('/data1/nkumar/RIPX/M_Files/')) %BL edit 

clear all
load RIPX_bath_guide
keyboard
load('RIPX_max_curlF_psi.mat','x_curlF','i_curlF','curlF');
%dirname        = '/data/ssuanda2/RIPX/QC1/';
dirname        = '/data1/ssuanda/RIPV/QC1/'; %BL modified 
dirname3       = '/data1/nkumar/RIPX/M_Files/RIPX_mstd_psi_curlF/';
fname          = 'qaqc_';
NS             = 9; % BL modified, original is 1 
NT             = 120;
rho0           = 1;
fac            = 1/16;
g              = 9.81;
sname          = 'RIPX_QC_';
ss_index       = 5;

xb                   = -QCstats(:,8);
Gamma                = QCstats(:,9)./QCstats(:,7);
Tp                   = 1./QCstats(:,4);
Linf                 = g*(Tp.^2)./(2*pi);
Lb                   = 2*pi./get_wavenum(2*pi./AA.A(:,7),QCstats(:,7));
uss                  = sqrt(g.*QCstats(:,7));
slp                  = AA.A(:,1);
stb                  = QCstats(:,3)./Lb;
stpin                = QCstats(:,3)./Linf;
Irbo                 = slp./sqrt(stpin);


%Define a structure:
S                    = struct; 

for i=NS:1:NT
    str        = sprintf('%02d',i);
    Fname      = [dirname,fname,str,'.mat'];
    eval(['load ','''',Fname,'''']);
    
    if i<10
        str      = sprintf('run_%d/',i);
        dirname2 = '/data1/ssuanda/RIPX/';
        fname2   = sprintf('mv_snap_%04d_eddy.mat',i);
        Fname2   = [dirname2,str,fname2];
        eval(['load ','''',Fname2,'''']);
        
    else
        str      = sprintf('run_%d/',i);
        dirname2 = '/data/ssuanda2/RIPX/';
        fname2   = sprintf('mv_snap_%04d_eddy.mat',i);
        Fname2   = [dirname2,str,fname2];
        eval(['load ','''',Fname2,'''']);
    end
    if i<115
        xpos = 270;
    else
        xpos = 320;
    end
    %std_visc = nanmean(sqrt(Bnvar(xpos:end,:)),2); %original
    std_visc = mean(sqrt(Bnvar(xpos:end,:)),2,"omitnan");% BL modified
    std_visc = std_visc(1:ss_index:end);
    
    %Get x and bathy again
    str0         = sprintf('%04d',i);
    eval(sprintf('bt = runlist.bath.bath_%s;',num2str(i)));
    eval(sprintf('X = bath.b0%s.x;',num2str(bt(9))));
    X            = X(twocol(i,1):twocol(i,2)); 
    X            = X(:);
    
    eval(sprintf('h = bath.b0%s.h;',num2str(bt(9))));
    h                = h(twocol(i,1):twocol(i,2)); 
    h                = h(:);
    
    % Get curlF_br
    fname3           = [dirname3,'RIPX_mstd_psi_curlF_',str0,'.mat'];
    eval(['load ',fname3]);
    curlF_mstd       = curlF_mstd(:);
    
    %Get Hs
    Hs               = nanmean(Hsig,2);
    Hs               = Hs(:);
    
    %Get dir spread
    sigma_th         = nanmean(dirsp,2);
    sigma_th         = sigma_th(:);
    
    %Get mean direction
    mean_th         = nanmean(mdir,2);
    mean_th         = mean_th(:);
    
    %Get a1 
    a1               = ((pi/180*sigma_th).^2/2)-1;
    
    %Get energy flux
    ECG              = fac*g*Hs.^2.*sqrt(g*H2(:)).*a1;
    dx               = 5;
    [~,dECG]         = calculate_gradient_from_central_difference(ECG,3,dx);
    dECG(end+1)      = 0;
    dECG             = [0 dECG];
    j                = find(X>0);
    dECG(j)          = 0;
    
    %Save variables
    S(i).X       = X;
    S(i).X2      = X2(:);
    S(i).h       = h;
    S(i).h2      = H2(:);
    S(i).curlF     = curlF(i);
    S(i).x_curlF   = x_curlF(i);
    S(i).i_curlF   = i_curlF(i);
    S(i).curlF_std = curlF_mstd;
    S(i).std_visc  = std_visc;
    S(i).xb        = xb(i);
    S(i).gamma_c   = Gamma(i);
    S(i).gamma_v   = Hs(:)./H2(:);
    S(i).Linf      = Linf(i);
    S(i).Lb        = Lb(i);
    S(i).m         = slp(i);
    S(i).stpb      = stb(i);
    S(i).stpin     = stpin(i);
    S(i).Irbo      = Irbo(i);
    S(i).Tp_T      = Tp(i);
    S(i).Hs        = Hs;
    S(i).sigma_th  = sigma_th;
    S(i).mean_th   = mean_th;
    S(i).a1        = a1;
    S(i).Hs_T      = QCstats(i,3);
    S(i).fp_T      = QCstats(i,4);
    S(i).sigma_T   = QCstats(i,6);
    S(i).hb        = QCstats(i,7);
    S(i).Hs_b      = QCstats(i,9);
    S(i).sigma_b   = QCstats(i,11);
    S(i).ECG       = ECG;
    S(i).dECG      = dECG;
    S(i).Hs_e      = AA.A(i,6);
    S(i).Tp_e      = AA.A(i,7);
    S(i).kh_e      = AA.A(i,8);
    S(i).sigma_e   = AA.A(i,10);    
    disp(i)
end

save qa_qc_RIPX_NK.mat S

figure
subplot(411)
plot(S(1).X2,S(1).Hs,'b','linewidth',3)
hold on
xlim([-300 0])
grid
plot(S(4).X2,S(4).Hs,'r','linewidth',3)

subplot(412)
plot(S(1).X,S(1).curlF_std,'b','linewidth',3)
hold on
xlim([-300 0])
grid
plot(S(4).X,S(4).curlF_std,'r','linewidth',3)

subplot(413)
plot(S(1).X2,S(1).std_visc,'b','linewidth',3)
xlim([-300 0])
grid
hold on
plot([-300:0],[-300:0]*0+0.01,'m');
plot(S(4).X2,S(4).std_visc,'r','linewidth',3)

subplot(414)
plot(S(1).X2,S(1).dECG,'b','linewidth',0.5)
xlim([-300 0])
dECG    = S(1).dECG;
j       = find(S(1).std_visc<0.01);
dECG(j) = NaN;
hold on
plot(S(1).X2,dECG,'b','linewidth',3)
grid
plot(S(4).X2,S(4).dECG,'r','linewidth',0.5)
dECG    = S(4).dECG;
j       = find(S(4).std_visc<0.01);
dECG(j) = NaN;
plot(S(4).X2,dECG,'r','linewidth',3)
 