clear all
close all

rnum_tot = 120;

CXT_ALL = cell(rnum_tot,1);

yindex = [100 200 300 400 500 600 700 800 900 1000 1100 1200];
ny2 = length(yindex);


nt_lag = 50;
ntlag = 2*nt_lag+1;
nx_lag = 10;

for rnum = 1:rnum_tot,
    disp(sprintf('starting rnum = %d',rnum));
    run = rnum;
  
    fname=sprintf('/data1/nkumar/RIPX/M_Files/RIPX_psi_curlF/RIPX_psi_curlF_%04d',run);

    load(fname,'curlF')

[nx,ny,nt]=size(curlF);

if (ny~=1200),
    disp('** ERROR. ny <> 1200. BAD!  abort');
    return;
end;

disp(sprintf('nt = %d',nt))
if (nt<4400)
    nt
    disp('** ERROR. nt <> 4801. BAD!  abort');
    return;
end;

CXT = zeros(nx,nx_lag,ntlag);
CXT_var = zeros(nx,nx_lag,ntlag);
CXT_G1_G1= zeros(nx,nx_lag,ntlag);
VG = zeros(nx,1);  % variance of G
cxt_tot = zeros(ntlag,1);
c_G2_G2_tot = zeros(ntlag,1);
cxt_tot_mag=0;
disp('starting looping')
for i=1:nx-nx_lag % 1--- xdim-10 --------> crosshore
    for j=1:nx_lag % 1--- 10 -------> x-lag %%%%%%%%%%%%%%%%%%%%%%%%%%% xlag
        cxt_tot = 0*cxt_tot;
        for k=1:ny2 %  --------> yavg
          iy = yindex(k);
          G1 = squeeze(curlF(i,iy,:)); %time series at given (x,y)
          G2 = squeeze(curlF(i+(j-1),iy,:)); %time series at given (x+xlag,y)
          G1 = G1-mean(G1); % ------------> curlF(x1,y1,t)
          G2 = G2-mean(G2); %--------------> curlF(x1+dx,y1,t+dt)
          cxt = xcorr(G2,G1,nt_lag); %BL modified
%           plot(G1)
%           hold on 
%           plot(G2)
%           hold off
%           drawnow
%           pause(0.5)
          %c_G1_G1 = xcorr(G1,G1,nt_lag);
          c_G2_G2 = xcorr(G2,G2,nt_lag); %  101 * 1 
          %cxt_mag = sqrt(c_G1_G1(51)*c_G2_G2(51)); % scalar
          %%         %%%%%%%%%%%%%%%%% for debug
%         if i == 130 && j == 8
%             keyboard;
%         end
        %%
          cxt_tot = sum(cat(3,cxt_tot,cxt),3,'omitnan');%BL modified
          c_G2_G2_tot = sum(cat(3,c_G2_G2_tot,c_G2_G2),3,'omitnan'); %BL modified

        end;
        cxt_tot = cxt_tot/ny2; %alongshore mean of covariance
        c_G2_G2_tot = c_G2_G2_tot/ny2;
        CXT_var(i,j,:) = cxt_tot;
        CXT_G1_G1(i,j,:) = c_G2_G2_tot; 
    end; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% x lag
    
    CXT(i,:,:) = CXT_var(i,:,:)./sqrt(CXT_var(i,1,51)*CXT_G1_G1(i,:,:));
end;

clear curlFbr;

CXT_ALL{rnum} = CXT;
CXT_ALL(rnum)
disp(sprintf('finished w rnum = %d',rnum));
end;
README = ['xcorr(G2,G1) and normalized using Falk sum convention'];

%save('/data1/bliu/data/raw/CXT_ALL_sumnormalized_BLrun.mat', 'CXT_ALL','README')
