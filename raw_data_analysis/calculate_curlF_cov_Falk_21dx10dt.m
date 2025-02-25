clear all
close all

rnum_tot = 120;

CXT_ALL = cell(rnum_tot,1);
CXT_ALL_var = cell(rnum_tot,1);


yindex = [100 200 300 400 500 600 700 800 900 1000 1100 1200];
ny2 = length(yindex);


nt_lag = 10;
ntlag = 2*nt_lag+1;
nx_lag = 21;

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
for i=1:nx-nx_lag % 1--- xdim-15 --------> crosshore
    for j=1:nx_lag % 1--- 15 -------> x-lag %%%%%%%%%%%%%%%%%%%%%%%%%%% xlag
        cxt_tot = 0*cxt_tot;
        for k=1:ny2 %  --------> yavg
          iy = yindex(k);
          G1 = squeeze(curlF(i,iy,:)); %time series at given (x,y)
          G2 = squeeze(curlF(i+(j-1),iy,:)); %time series at given (x+xlag,y)
          G1 = G1-mean(G1); % ------------> curlF(x1,y1,t)
          G2 = G2-mean(G2); %--------------> curlF(x1+dx,y1,t+dt)
          cxt = xcorr(G2,G1,nt_lag); %BL modified
          cxt_tot = sum(cat(3,cxt_tot,cxt),3,'omitnan');%BL modified: add cxt for diff y loc for avg

        end
        cxt_tot = cxt_tot/ny2; %alongshore mean of covariance
        CXT_var(i,j,:) = cxt_tot;
    end; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% x lag
%              %%         %%%%%%%%%%%%%%%%% for debug
% 
%         test = squeeze(CXT_var(i,:,:));
%         [test_max,max_ind] = max(test,[],"all","linear");
%         [row,col] = ind2sub(size(test),max_ind);
%         report = (row,col)
%         if ((row~=1)|(col ~=51)) && abs(test_max)>0.01
%             keyboard
%         end 
% 
%         
%         %%
    c0 = CXT_var(i,1,nt_lag+1);
    CXT(i,:,:) = CXT_var(i,:,:)/c0;
end;

clear curlFbr;

CXT_ALL{rnum} = CXT;
CXT_ALL_var{rnum} = CXT_var;
disp(sprintf('finished w rnum = %d',rnum));
end;
README = ['xcorr(G2,G1) and normalized using Falk sum convention and also include variance and this increased dx domain and decreased dt domain'];

save('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt.mat', 'CXT_ALL','CXT_ALL_var','README')
