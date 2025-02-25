clear all
addpath(genpath('/home/n2kumar/MAT_Mytoolbox/jlab/'));
NT   = 120;

for i=2:1:NT
    str0     = sprintf('%04d',i);
    fname1   = ['RIPX_psi_curlF_',str0,'.mat'];
    eval(['Fname1=matfile(','''',fname1,'''',');']);
    
    [p,nc,r] = size(Fname1.curlF); % (x, y, t)
    ntapy    = 5;
    ppy      = (ntapy+1)/2;
    [psic,~] = sleptap(nc,ppy);
    nk       = 2;
    for j = 1:1:p %cross-shore index
    %for j = 120:1:p %BL modified for testing (need to change otherwise will lead to bug)

        count    = 0; %time index
        for k = 1:nk:r %time index
            count                     = count+1;
            cm                        = nanmean(squeeze(Fname1.curlF(j,:,k))); % mean in y
            c                         = squeeze(Fname1.curlF(j,:,k))-cm; %demeaned y slice
            [fck,Sck(j,:,count)]      = mspec(c',psic);
            if mod(k,101)==0
                disp(k)
            end
        end
        disp(j)
    end
    disp('Done with alongshore wavenumber part')
    clear j k c psic
    disp(['At file number',str0])
    sname = ['RIPX_curlF_spec_',str0,'.mat'];
    eval(['save -v7.3 ',sname,' fck',' Sck'])
end

    
%     [p,q,nt] = size(Fname1.curlF);
%     ntapy    = 5;
%     ppy      = (ntapy+1)/2;
%     [psit,~] = sleptap(nt,ppy);
%     count    = 0;
%     for j = 1:dx:p
%         count = count+1;
%         count2= 0;
%         for k = 1:dy:q
%             count2           = count2+1;
%             cm               = nanmean(squeeze(Fname1.curlF(j,k,:)));
%             c                = squeeze(Fname1.curlF(j,k,:))-cm; 
%             [fct,Sct(count,count2,:)] = mspec(c,psit);
%         end
%         disp(j)
%     end    
%     
    
%     
% [~,r,nu]= size(U);
% [~,~,nv]= size(V);
% [~,~,nw]= size(W);
% [~,~,nt]= size(T);
% 
% ntapy    = 5;
% ppy      = (ntapy+1)/2;
% [psiu,~] = sleptap(nu,ppy);
% [psiv,~] = sleptap(nv,ppy);
% [psiw,~] = sleptap(nw,ppy);
% [psit,~] = sleptap(nt,ppy);
% 
% for i=1:1:p
%     count = 0;
%     for j=1:1:r
%         um = nanmean(squeeze(U(i,j,:)));
%         if count>0
%             vm = nanmean(squeeze(V(i,count,:)));
%         end
%         wm = nanmean(squeeze(W(i,j,:)));
%         tm = nanmean(squeeze(T(i,j,:)));
%         
%         u  = squeeze(U(i,j,:))-um;
%         [fut,Sut(i,j,:)]=mspec(u,psiu);
%         
%         if count>0
%             v                    = squeeze(V(i,count,:))-vm;
%             [fvt,Svt(i,count,:)] = mspec(v,psiv);
%         end
%         
%         w  = squeeze(W(i,j,:))-wm;
%         [fwt,Swt(i,j,:)]=mspec(w,psiw);
%         
%         t  = squeeze(T(i,j,:))-tm;
%         [ftt,Stt(i,j,:)]=mspec(t,psit);
%         count = count+1;
%     end
%     disp(i)
% end

