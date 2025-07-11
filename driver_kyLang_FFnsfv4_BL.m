% Bingchen Liu May 27, 2025
% modified from Matt's original code 'driver_kyLang_FFnsfv4'


addpath(genpath('/home/n2kumar/MAT_Mytoolbox/jlab/'));

clear all
close all

%%
Ly=1500;
dy=1;
Lx=350;
dx=1;
x=0:dx:Lx;
nx=length(x);
y=0:dy:Ly;
ny=length(y);

%%
[X,Y]=meshgrid(x,y);
%%
if mod(ny,2)==0
    ky=[-ny/2+1:1:ny/2 ]';    
else
    ky=[-(ny-1)/2:1:(ny-1)/2 ]';
end
ky=ky/Ly;
dky=ky(2)-ky(1);

%% get Weibull spectra
ky0=.02;%peak alongshore wavenumber
beta= 1.4; %modified gamma in the new paper version 
%lam = ky0*beta/(beta-1)^(1/beta);
lam = ky0/((beta-1)/beta)^(1/beta); %BL: I think this is the correction version
Aky = (beta/(lam^beta))*abs(ky).^(beta-1).* exp(-(abs(ky/lam)).^(beta)); % Weibull spectra
var0=sum(Aky)*dky; % the variance I think I should get
Aky = Aky / var0; % the variance should now be 1
% becasue int (Aky dky) = 1 here

phi=2*pi*rand(size(Aky)); %random phase 
Ak=(2*ny*Aky).^.5.*exp(1i*phi); % Fourier component with random phase; Ny is for Matlab fft scaling

%% generate spatial series using the Weibull spectra 
clf
psi0=real(ifft(ifftshift(Ak))); % inverse fft back to space domain and take the real part
varpsi=var(psi0); % the total variance I actually get, this is now 1

figure(1)
subplot(3,1,1)
plot(y,psi0)
%[ky2,Syy]=fft_data(psi0,dy);
[psic,~] = sleptap(length(psi0),3); % 3 is related to number of taper --> ntaper = 3*2-1
[ky2,Syy] = mspec(psi0,psic); %output of ky2 is in 2*pi/T or wavelength
ky2 = ky2/(2*pi); %convert to cpm
subplot(3,1,2)
plot(ky,Aky,ky2,Syy)
subplot(3,1,3)
loglog(ky,2*Aky,ky2,Syy)


%% AA: Fourier component (make 2D via AR1 in x) 


nky=size(Aky,1);
AA= zeros(nky,nx); % 2D field 
AA(:,1)=Ak;
LX=5;
tau = 0.75;%decorrelation time scale
a1= exp(-1/tau);%a1 in discrete AR1, how much info to take from previous step 
%LX= 1; % BL test 
c=.5*(4*LX^2-2*dx*LX+dx^2)/(2*LX^2-dx*LX); %scalar 
cc=(2*c*(2*ny*Aky)*dx/LX).^.5 ; %loglog(ky,Aky) is the spectra 

for b=2:length(x)
    AA(:,b)=AA(:,b-1)-AA(:,b-1)/LX * dx + cc * 2^(-.5) .* (randn(nky,1) + (-1)^.5*randn(nky,1)); % add white noise that has Weibull spectra shape (information in cc) --matt version
    %AA(:,b)=a1*AA(:,b-1) + cc * 2^(-.5) .* (randn(nky,1) + (-1)^.5*randn(nky,1)); %BL version (only works for dx=1)
    %AA(:,b)=AA(:,b-1)-AA(:,b-1)/LX * dx + cc * 2^(-.5) .* (randn(nky,1)); %BL testing 

end

% AA is the 2D field now with (:,i) being the alongshore 
%%
% nt=1e5;
% X=zeros(1,nt);
% vx=.01;
% tau=10;
% X(1)=vx.^.5 * randn(1);
% dt=1;
% for b=2:nt
%     X(b)=X(b-1)-X(b-1)*dt/tau + ( 2 * vx * dt / tau ).^.5 * randn(1);
% end
%% BL test 
cc=(2*c*(2*ny*Aky)*dx/LX).^.5 ;
noise_Fcomp = cc * 2^(-.5) .* (randn(nky,1) + (-1)^.5*randn(nky,1));
noise_spac  = real(ifft(ifftshift(noise_Fcomp)));
[f_noise,S_noise] = fft_data(noise_spac,dy);
%% AA--> Psi0 (convert to physical space )
for b=1:length(x)
    Psi0(:,b)=real(ifft(ifftshift(AA(:,b)))); 
    %Psi0(:,b)=abs(ifft(ifftshift(AA(:,b))));  %BL test

end

[f1,S1] = fft_data(Psi0(:,1),dy);
[f2,S2] = fft_data(Psi0(:,330),dy);
figure()
loglog(f1,S1,f2,S2,ky,2*Aky,f_noise,S_noise)
legend('loc1','loc2','Theory','noise')
%%
mxx0=.075;
x00=25;
bb=x00/2^(1/3);
aa=mxx0*2^(1/3)*1.5/bb;
Psim1=aa*X./(1+(X/bb).^3); %G0 magnitude function 
x01=30;
x02=130;
lxx=1.5;
Psim2=0*Psim1; %make Width function pre-allocated with 0
wlength=60;
x000=[30:wlength:350];

x000b=[40:5:max(x)-40];

W02_1 = wlength / pi^.5 / lxx;
dum = exp(.25*lxx^2/LX^2) * ( erf( (lxx^2+wlength*LX)/(2*lxx*LX)) - erf(lxx/(2*LX)) );
W02_2i = dum;
W02 = W02_1 / W02_2i;

ct =3; %estimate for c\tau
%A1 = 8.39;
%A2 = 

for a=1:length(x000b) %generate multiple width functions (more than needed to validate) in x individually (not combineed) --> check normalization for all of the humps
    Psim2=(W02_1*exp(-(X-x000b(a)).^2 / lxx^2 ) ).^.5; %matt format
    %Psim2=(A1/cos(A4))*exp(-A2*abs(X-x000b(a)))*cos(A3*abs(X-x000b(a))+A4); %BL version    
    %Psim2=(W02*exp(-(X-x000b(a)).^2 / lxx^2 ) ).^.5;
    ii=find(x>=x000b(a)-.5*wlength & x<=x000b(a)+.5*wlength); %get info around 
    dum=Psim2(:,ii).*Psi0(:,ii);
    msqdum(a)=mean(dum.^2,'all'); % this should equal 1 for all a same as var(dum)=1
    msqdum0(a)=mean(Psim2(:,ii).^2,'all'); %should =1, same as var(Psi0)=1
end

figure(2)
clf
plot(msqdum,'o')
hold on
plot(1:a,mean(msqdum)*ones(1,a),'--k')
plot(msqdum0,'x')
xlabel('strip number')
ylabel('$\langle [W(x-ct) \psi(x-ct,y)]^2 \rangle$','Interpreter','latex')
legend('individual strips','mean','just \psi')
title('Normaliztion Check')


Psim2=0*Psim2;
for a=1:length(x000) %combine all the humps 
    Psim2=Psim2+(wlength/pi^.5/lxx*exp(-(X-x000(a)).^2 / lxx^2 ) ).^.5;
    ii=find(x>=x000(a)-.5*wlength & x<=x000(a)+.5*wlength);
end
Psi=Psim1.*Psim2.*Psi0; %final field. 
% Psi0: background forcing with given Swb
% Psim1: magnitude from G0
% Psim2: width/envolope function 

%%
%paperfigure;
figure()
subplot(1,2,1)
contour(x,y,Psi,24)
axis equal
axis xy
colorbar
subplot(1,2,2)
contour(x,y,Psim1.*Psim2,24)
axis equal
axis xy
colorbar

%%
%paperfigure;
figure()
subplot(2,2,1)
pcolor(x,y,Psi0)
mxx2=max(max(abs(Psi0)));
mxx=.95*mxx2;
caxis(mxx*[-1 1])
colormap(bluewhitered)
shading flat
axis xy
colorbar
set(gca,'xlim',[0 400],'ylim',[0 1.5e3])
subplot(2,2,2)
pcolor(x,y,Psim1)
caxis([0 mxx0])
colormap(bluewhitered)
shading flat
axis xy
colorbar
%dum=Psi2;
%dum(abs(Psi2)<.4*mxx0)=0;
set(gca,'xlim',[0 400],'ylim',[0 1.5e3])
subplot(2,2,3)
pcolor(x,y,Psim2)
caxis([0 1])
colormap(bluewhitered)
shading flat
axis xy
colorbar
set(gca,'xlim',[0 400],'ylim',[0 1.5e3])
subplot(2,2,4)
pcolor(x,y,Psi)
caxis([-.1 .1])
colormap(bluewhitered)
shading flat
axis xy
colorbar
set(gca,'xlim',[0 400],'ylim',[0 1.5e3])
title('$G = \nabla \times F_{\mathrm{br}}$')

%% proposal fig <=========================
%paperfigure;
x00=.15;
y00=.35;
ly2=.225;
dy1=.025;
dx1=.08;
lx1=(1-2*x00-dx1)/2;
ly1=(ly2-dy1)/2;
axp(2,:)=[x00 y00 lx1 ly1];
axp(1,:)=axp(2,:)+[0 ly1+dy1 0 0];
axp(3,:)=axp(2,:)+[lx1+dx1 0 0 -ly1+ly2];
figure
axh(1)=axes('position',axp(1,:));
plot(-x,Psim1(1,:),'k')
set(gca,'xlim',[-200 0],'ylim',[0 .08],'layer','top','xticklabel',[])
ylabel('$G_0$ [s$^{-2}$]')
rxl=.03;
ryl=.85;
%txth=label_panel('(a)',rxl,ryl);
axh(2)=axes('position',axp(2,:));
plot(-x,Psim2(1,:),'k')
set(gca,'xlim',[-200 0],'ylim',[0 (wlength/pi^.5/lxx)^.5*1.1],'layer','top')
ylabel('$W$')
xlabel('$x$ [m]')
%txth=label_panel('(b)',rxl,ryl);
axh(3)=axes('position',axp(3,:));
pcolor(-x,y,Psi)
caxis([-.25 .25])
colormap(bluewhitered)
colorbar
shading interp
ylabel('$y$ [m]')
xlabel('$x$ [m]')
title('$G = \nabla \times F_{\mathrm{br}}$')
set(gca,'xlim',[-200 0],'ylim',[0 250],'layer','top')
grid on
%txth=label_panel('(c)',rxl,.92);

return


%%
x=linspace(-10,10,1001)';
lx=1.5;
LX=5;
d = [-5:2.5:5];

[D,X]=meshgrid(d,x);

W1=(1 / pi^.5 / lx * exp( - (D/lx).^2 ) ).^.5;
W2=(1 / pi^.5 / lx * exp( - ((X+D)/lx).^2) ).^.5;
C =(exp(-abs(X)/LX));

WWC = W1.*W2.*C;

%%
figure(1)
clf
hold on
clrs=jet(length(d));
clrs=1-clrs;
for a=1:length(d)
    plot(x,WWC(:,a),'color',clrs(a,:))
end
plot(x,WWC(:,3),'--k','linewidth',4)
plot(x,max(WWC(:,3))*C(:,1),'g','linewidth',5)
legend('\phi=-5','-2.5','0','2.5','5','0','C')
ylabel('$W(\phi)W(\phi+\Delta)C(\Delta)$')
xlabel('$\Delta$')
title(['$W^2(\phi) = \pi^{-1/2} l_x^{-1} \exp(-\phi^2/l_x^2)$, $C=\exp(-| \Delta | / \tau_x )$, $l_x=$' num2str(lx) ', $\tau_x=$' num2str(LX)])

%%
figure(4)
set(gcf,'paperpositionmode','auto','color','w','renderer','painters')
export_fig fig_2021_FFnsf_curlF.pdf -pdf
!mv fig_2021_FFnsf_curlF.pdf /Users/mspydell/Desktop/

%%
figure(2)
set(gcf,'paperpositionmode','auto','color','w','renderer','painters')
export_fig fig_2021_FFnsf_varWpsi.pdf -pdf
!mv fig_2021_FFnsf_varWpsi.pdf /Users/mspydell/Desktop/
    

%%







return
