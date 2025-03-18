% Bingchen Liu Nov 15,2024
% Include different normalization using Rxy = <x y>/ {std(x)std(y)}
%
% modified on Jan 8, 2025 to include uncertainty range for CXT 

clear all
close all

rnum_tot = 120;
z= 1.96; %parameter for confidence interval 1.96 -- 95%

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
%CXT_var = zeros(nx,nx_lag,ntlag);
CXT_lowerbd= zeros(nx,nx_lag,ntlag);
CXT_upperbd= zeros(nx,nx_lag,ntlag);
CXT_uncert= zeros(nx,nx_lag,ntlag);
CXT_std= zeros(nx,nx_lag,ntlag);


CXT_G1_G1= zeros(nx,nx_lag,ntlag);
VG = zeros(nx,1);  % variance of G
cxt_tot = zeros(ntlag,1);
c_G2_G2_tot = zeros(ntlag,1);
cxt_tot_mag=0;
disp('starting looping')
for i=1:nx %-nx_lag % 1--- xdim-15 --------> crosshore ---> BL modified Feb 13,2025
    for j=1:nx_lag % 1--- 15 -------> x-lag %%%%%%%%%%%%%%%%%%%%%%%%%%% xlag
        cxt_tot = 0*cxt_tot;
        norm_cxt_tot=0;
        for k=1:ny2 %  --------> yavg
          
          
            % if rnum==2&&i==1&&k==1
            %     keyboard
            % end

          iy = yindex(k);
          G1 = squeeze(curlF(i,iy,:)); %time series at given (x,y)
          if i+j-1<=nx %BL modified Feb 18, 2025
              G2 = squeeze(curlF(i+(j-1),iy,:)); %time series at given (x+xlag,y)
          else
              G2 = nan(4801,1); %BL modified 
          end
          
          G1 = G1-mean(G1,'omitmissing'); % ------------> curlF(x1,y1,t)
          G2 = G2-mean(G2,'omitmissing'); %--------------> curlF(x1+dx,y1,t+dt)
          [cxt,~] = xcorr(G2,G1,nt_lag); %BL modified
          r_11=xcorr(G1,G1,nt_lag);
          r_22=xcorr(G2,G2,nt_lag);
          norm_cxt = sqrt(r_11(11)*r_22(11));
          cxt_tot = sum(cat(3,cxt_tot,cxt),3,'omitnan');%BL modified: add cxt for diff y loc for avg
          norm_cxt_tot=sum(cat(1,norm_cxt_tot,norm_cxt),1,'omitnan');
          cxt_aty(:,k) =cxt./norm_cxt; 
        end  %  --------> yavg
        cxt_tot = cxt_tot/ny2; %alongshore mean of covariance
        norm_cxt_tot = norm_cxt_tot/ny2;
        CXT_temp = cxt_tot/norm_cxt_tot;
        CXT(i,j,:) = CXT_temp;
        %%%%%%%%%%%%%%%% include upper and lower bound uncertainty estimate
        %%%%%%%%%%%%%%%% for CXT
% if i==120&&j==21
%     keyboard
% end 
        CXT_lowerbd(i,j,:) = ((1+CXT_temp)-(1-CXT_temp).*exp(2*z/sqrt(length(G1)*ny2/3-3)))./((1+CXT_temp)+(1-CXT_temp).*exp(2*z/sqrt(length(G1)*ny2/3-3))); %use z_0.25 = 1.96 --> 95% confidence interval
        CXT_upperbd(i,j,:) = ((1+CXT_temp)-(1-CXT_temp).*exp(-2*z/sqrt(length(G1)*ny2/3-3)))./((1+CXT_temp)+(1-CXT_temp).*exp(-2*z/sqrt(length(G1)*ny2/3-3))); %use z_0.25 = 1.96 --> 95% confidence interval
        CXT_uncert(i,j,:) =CXT_upperbd(i,j,:)-CXT_lowerbd(i,j,:); 
        CXT_std(i,j,:) = std(cxt_aty,0,2,'omitmissing');
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
    %c0 = CXT_var(i,1,nt_lag+1);
    %CXT(i,:,:) = CXT(i,:,:);
end;

clear curlF;

CXT_ALL{rnum} = CXT;
CXT_lowbd_ALL{rnum} = CXT_lowerbd;
CXT_upbd_ALL{rnum} = CXT_upperbd;
CXT_uncert_ALL{rnum} = CXT_uncert;
CXT_std_ALL{rnum} = CXT_std;
%CXT_ALL_var{rnum} = CXT_var;
disp(sprintf('finished w rnum = %d',rnum));
end;
README = ['generated in calculate_curlF_cov_Bingchen_21dx10dt, xcorr(G2,G1) and normalized using Bingchen sum convention and also include variance and this increased dx domain (20m) and decreased dt domain (10s), this also include upper and lower bound for CXT'];

save('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen_includeshore.mat', 'CXT_ALL','CXT_lowbd_ALL','CXT_upbd_ALL','CXT_std_ALL','CXT_uncert_ALL','README')
%save('/data1/bliu/data/raw/CXT_ALL_norm_and_var_21dx10dt_Bingchen.mat', 'CXT_ALL','CXT_lowbd_ALL','CXT_upbd_ALL','CXT_std_ALL','CXT_uncert_ALL','README')
