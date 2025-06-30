% Bingchen Liu June 27,2025
% This code solves the equation in the width function calculation in RIPX

C1=3.57;
C2=1.82;
C3 = 1.24;
ct=3;
A2=C1/ct;
A3=C2/ct;

syms theta x
alpha = exp(-2*3.57*theta)*cos(1.82*theta+x)*cos(1.82);
beta = exp(-2*3.57*theta)*cos(1.82*theta+x)*sin(1.82);

eqn=x+atan2(int(beta,theta,0,inf),int(alpha,theta,0,inf))==1.24;
sol= solve(eqn,x)

B1 = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+sol)*cos(1.82),theta,0,inf);%alpha
B2 = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+sol)*sin(1.82),theta,0,inf); %beta

A4 = sol;

C3_validate = double(sol+atan2(B2,int1))
A1 = (cos(A4))/(cos(C3)^0.5*(B1^2+B2^2)^0.25);

%% plot W function
close all 
C1= 3.57;
C2 = 1.82;
ct = 3; 
A2 = 3.57/ct;
A3 = 1.82/ct;
%A4=-A4;


dx=0.01;
x=-10:dx:10;
W = double((A1/cos(A4)).*exp(-A2*abs(x)).*cos(A3*abs(x)+A4));
int_W=  sum(W.^2)*dx

W1= exp(-A2*abs(x));
W2= cos(A3*abs(x)+A4);

CX = (1/cos(C3)).*exp(-C1*abs(x/ct)).*cos(C2*abs(x/ct)+C3);

figure()
plot(x,W,'LineWidth',2)
hold on 
plot(x,CX,'LineWidth',2)
hold on 
plot(x,W1,'LineWidth',2)
hold on 
plot(x,W2,'LineWidth',2)
hold off 
legend('W(tot)','C_X','exp','oscillation')
xlabel('$\theta$','Interpreter','latex')
ylabel('W')


[acf,lags] = xcorr(W,'unbiased');

figure()
plot(lags*dx,acf)
hold on 
plot(x,CX)
hold off 

%% Use Falk method for solving the equation
clear
syms theta 
A4_solve= 0:0.005:5;
for ind= 1:length(A4_solve)
B1_solve = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+A4_solve(ind))*cos(1.82),theta,0,inf);%alpha
B2_solve = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+A4_solve(ind))*sin(1.82),theta,0,inf); %beta
C3_solve(ind) = A4_solve(ind)+atan2(B2_solve,B1_solve);
end 

C3_solve = double(C3_solve);
[val,ind] = min(abs(C3_solve-1.24));
A4 = A4_solve(ind);
figure()
plot(A4_solve,C3_solve)

%% plot Falk method solution
close all 
C1= 3.57;
C2 = 1.82;
ct = 3; 
A2 = 3.57/ct;
A3 = 1.82/ct;
B1 = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4)*cos(A3*theta),theta,0,inf);%alpha
B2 = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4)*sin(A3*theta),theta,0,inf); %beta
A1 = (cos(-A4))/(cos(C3)^0.5*(B1^2+B2^2)^0.25);


dx=0.01;
x=-10:dx:10;
W = double((A1/cos(A4)).*exp(-A2*abs(x)).*cos(A3*abs(x)+A4));
int_W=  sum(W.^2)*dx;

CX = (1/cos(C3)).*exp(-C1*abs(x/ct)).*cos(C2*abs(x/ct)+C3);

[acf,lags] = xcorr(W,'unbiased');

figure()
plot(lags*dx,acf)
hold on 
plot(x,CX)
hold off 







