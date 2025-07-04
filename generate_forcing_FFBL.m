% Bingchen Liu June 27,2025
% This code generate the vorticity forcing using Falk's method with
% cross-shore: uniform forcing  
% alongshore: Fourier components in alongshore
% Time: AR1 in time in 3rd dimension

% Forcing decomposition G = Psi0*Psim1*Psim2 
% Psi0: background forcing with given Swb
% G0: magnitude from G0
% W: width/envolope function


addpath(genpath('/home/n2kumar/MAT_Mytoolbox/jlab/'));

clearvars -except test
close all

%% initial setup
Ly=1500;
dy=1;
Lx=350;
dx=1;
t_tot= 80;
dt=1;
x=0:dx:Lx;
nx=length(x);
y=0:dy:Ly;
ny=length(y);
t=0:dt:t_tot;
nt = length(t);
h = 2; %water depth 
g= 9.81; %gravitational constant 
c_phase = -sqrt(g*h);



[X,Y]=meshgrid(x,y);

if mod(ny,2)==0
    ky=[-ny/2+1:1:ny/2 ]';    
else
    ky=[-(ny-1)/2:1:(ny-1)/2 ]';
end
ky=ky/Ly;
dky=ky(2)-ky(1);


%% get the magnitude function G0 first
mxx0=.075;
x00=120;
bb=x00/2^(1/3);
aa=mxx0*2^(1/3)*1.5/bb;
G0=aa*(X)./(1+((X)/bb).^3); %G0 magnitude function
G0 = repmat(G0,1,1,nt);
plot(x,G0(1,:,1))

%% get Weibull spectra
ky0=.02;%peF_comp alongshore wavenumber
beta= 1.4; %modified gamma in the new paper version 
%lam_matt = ky0*beta/(beta-1)^(1/beta);
lam = ky0/((beta-1)/beta)^(1/beta); %BL: I think this is the correction version
S_wb = (beta/(lam^beta))*abs(ky).^(beta-1).* exp(-(abs(ky/lam)).^(beta)); % Weibull spectra
%S_wb_matt = (beta/(lam_matt^beta))*abs(ky).^(beta-1).* exp(-(abs(ky/lam_matt)).^(beta)); % Weibull spectra

var0=sum(S_wb)*dky; % the variance I think I should get
S_wb = S_wb / var0; % the variance should now be 1
% becasue int (S_wb dky) = 1 here

phi=2*pi*rand(size(S_wb)); %random phase 
F_comp=(2*ny*S_wb).^.5.*exp(1i*phi); % Fourier component with random phase; Ny is for Matlab fft scaling

%% generate spatial series using the Weibull spectra 

psi0_test=real(ifft(ifftshift(F_comp))); % inverse fft back to space domain and tF_compe the real part
varpsi=var(psi0_test) % the total variance I actually get, this is now 1

% ===================> verify the spectra has the required Weibull shape
figure(1)
subplot(3,1,1)
plot(y,psi0_test)
%[ky2,Syy]=fft_data(psi0,dy);
[psic,~] = sleptap(length(psi0_test),3); % 3 is related to number of taper --> ntaper = 3*2-1
[ky2,Syy] = mspec(psi0_test,psic); %output of ky2 is in 2*pi/T or wavelength
ky2 = ky2/(2*pi); %convert to cpm
subplot(3,1,2)
plot(ky,S_wb,ky2,Syy)
subplot(3,1,3)
loglog(ky,2*S_wb,ky2,Syy)



%% AA: using F_comp: Fourier component w random phase (AR1 in time, fix normalization!!!!!!!) 


nky=size(S_wb,1);
AA= zeros(nky,nx,nt); % 2D field 
AA(:,1,1)=F_comp; %initilize the first alongshore forcing strip @t=0 using Swb Fourier component with random phase 

% ================> AR1 input <==========================
tau = 0.75;%decorrelation time scale
a1= exp(-1/tau);%a1 in discrete AR1, how much info to tF_compe from previous step 


% =================> normalization for AR1 <=====================
LX=5;
c=.5*(4*LX^2-2*dx*LX+dx^2)/(2*LX^2-dx*LX); %scalar 
cc=(2*c*(2*ny*S_wb)*dx/LX).^.5 ;  %S_wb^(0.5) but with different normalization


% =========> test noise 
noise_Fcomp = cc * 2^(-.5) .* (randn(nky,1) + (-1)^.5*randn(nky,1));
noise_spac  = real(ifft(ifftshift(noise_Fcomp)));

[ky_noise,S_noise] = fft_data(noise_spac,dy);
figure()
loglog(ky_noise,S_noise)
hold on 
loglog(ky,2*S_wb)
hold on 
loglog(ky,cc)
hold off 
legend('noise','S_{WB}','cc')


% ===================> AR1 in time <=======================
%create a slice in a given xloc and span in y and t
factor = sqrt(2/tau); %BL artificial factor to figure out the normalization.it's just a test. need to calculate the real one
for tind = 2:nt
    AA(:,1,tind)=a1*AA(:,1,tind-1)+ cc*factor * 2^(-.5) .* (randn(nky,1) + (-1)^.5*randn(nky,1)); %AR1 it with coef a1 and cc  
    % cc: 'noise' Fourier component that satisfies Weibull spectra,
    % different normalization based on a1 
end 

% =================> fill in AA with constant strip @x=1 at different time
%then copy the above slice to all xloc 
%note using uniform Fourier component in x 
for xind = 2:nx
    AA(:,xind,:)=AA(:,xind-1,:);
end 


%% AA--> Psi0(background forcing phase):convert to physical space 

for tind=1:nt
    for xind = 1:nx
        Psi0(:,xind,tind)=real(ifft(ifftshift(AA(:,xind,tind)))); %forcing in physical space 
    end 
end

%============================> visual test (forcing look and variance )
%%%%%%%%%%%%%%% forcing visual check
% for tind = 1:100
%     imagesc(squeeze(Psi0(:,:,tind)))
%     drawnow
%     pause(0.1)
% end 

%%%%%%%%%%%%%% forcing spectra check 
[f1,S1] = fft_data(Psi0(:,1,1),dy);
[f2,S2] = fft_data(Psi0(:,1,80),dy);

figure()
loglog(f1,S1,f2,S2,ky,2*S_wb)
legend('loc1','loc2','Theory')
xlabel('k_y(m^{-1})')
ylabel('Spectra energy')


%%%%%%%%%%%% forcing variance check (var(Psi0) ~1)

for i=1:nt
    var_Psi0(i) = var(squeeze(Psi0(:,1,i)));
end 
figure
plot(var_Psi0)




%% W (width function)

%==========================> width parameters
ct = 3; %estimate for c\tau
A1 = 1.82;
A2 = 3.57/ct;
A3 = 1.82/ct;
A4 = 1.1;
wlength =60; %crest separation 
x0_start = 0;
for tind=1:nt
    x0_start= floor(mod(x0_start+c_phase*dt,wlength));%location of start of the crest 
    x0{tind}=(wlength/2)+x0_start:wlength:Lx;
end 

x000b=[40:5:max(x)-40]; %location between 40~300 for normalization testing 
   

%%%%%%%%%%%%%%%%%%%%%%%% normalization check 
% for a=1:length(x000b) %generate multiple width functions (more than needed to validate) in x individually (not combineed) --> check normalization for all of the humps
%     %Psim2=(W02_1*exp(-(X-x000b(a)).^2 / lxx^2 ) ).^.5; %matt format
%     Psim2=(A1/cos(A4)).*exp(-A2.*abs(X-x000b(a))).*cos(A3.*abs(X-x000b(a))+A4); %BL version    
% 
%     % ii=find(x>=x000b(a)-.5*wlength & x<=x000b(a)+.5*wlength); %get info around 
%     % dum=Psim2(:,ii).*Psi0(:,ii);
%     % msqdum(a)=mean(dum.^2,'all'); % this should equal 1 for all a
%     % msqdum0(a)=mean(Psi0(:,ii).^2,'all');
% end

%==========================> generate multiple width function in the
%xdomain  
W=0*Psi0;
for tind = 1:nt
for a=1:length(x0{tind}) %combine all the humps 
    W(:,:,tind)=W(:,:,tind)+(A1/cos(A4)).*exp(-A2.*abs(X-x0{tind}(a))).*cos(A3.*abs(X-x0{tind}(a))+A4);
    %ii=find(x>=x000(a)-.5*wlength & x<=x000(a)+.5*wlength);
end
end 
cFbr=G0.*W.*Psi0; %final field.

head = 'data generated from generate_forcing_FFBL';

save('/data1/bliu/data/parameterization_example','cFbr','G0','W','X','Y',"head")
%% %%%%%%%%%%%%%%%%%%% visual of the width function 

figure(10)
for tind = 1:nt
    subplot(221)
    plot(-X(1,:),G0(1,:,tind))
    title('G0')
    subplot(222)
    plot(-X(1,:),W(1,:,tind))
    title('W')
    subplot(223)
    pcolorcen(-X,Y,Psi0(:,:,tind))
    cmocean('balance');
    caxis([-5,5])
    title('$\tilde{a}*b$','Interpreter','latex')
    subplot(224)
    pcolorcen(-X,Y,cFbr(:,:,tind))
    cmocean('balance');
    caxis([-0.2,0.2])
    title('$\nabla \times F_{br}$','Interpreter','latex')
    drawnow
    pause(0.5)
end 

%%
figure()
pcolorcen(-X,Y,cFbr(:,:,10))
cmocean('balance');
caxis([-0.3,0.3])

