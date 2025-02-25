clear all

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
  
    fname=sprintf('../RIPX_psi_curlF/RIPX_psi_curlF_%04d',run);

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
VG = zeros(nx,1);  % variance of G
cxt_tot = zeros(ntlag,1);
disp('starting looping')
for i=1:nx-nx_lag,
    for j=1:nx_lag,
        cxt_tot = 0*cxt_tot;
        for k=1:ny2,
          iy = yindex(k);
          G1 = squeeze(curlF(i,iy,:));
          G2 = squeeze(curlF(i+(j-1),iy,:));
          G1 = G1-mean(G1);
          G2 = G2-mean(G2);
          cxt = xcorr(G1,G2,nt_lag);
          cxt_tot = cxt_tot+cxt;
        end;
        cxt_tot = cxt_tot/ny2;
        CXT(i,j,:) = cxt_tot;
    end;
end;

clear curlFbr;

CXT_ALL{rnum} = CXT;
CXT_ALL(rnum)
disp(sprintf('finished w rnum = %d',rnum));
end;

save CXT_ALL.mat CXT_ALL
