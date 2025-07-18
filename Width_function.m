% Bingchen Liu June 30,2025
% This code solves the equation in the width function calculation in RIPX
% added Falk's method in solve for A4

C1=3.57;
C2=1.82;
C3 = 1.24;
ct=1;
A2=C1/ct;
A3=C2/ct;

syms theta x
B1_integrand = exp(-2*A2*theta)*cos(A3*theta+x)*cos(A3*theta);
B2_integrand = exp(-2*A2*theta)*cos(A3*theta+x)*sin(A3*theta);

eqn=x+atan2(int(B2_integrand,theta,0,inf),int(B1_integrand,theta,0,inf))==1.24;
sol= solve(eqn,x)

B1 = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+sol)*cos(1.82*theta),theta,0,inf);%alpha
B2 = 2*int(exp(-2*3.57*theta)*cos(1.82*theta+sol)*sin(1.82*theta),theta,0,inf); %beta

A4 = sol;

C3_validate = double(sol+atan2(B2,int1))
A1 = (cos(A4))/(cos(C3)^0.5*(B1^2+B2^2)^0.25);

% plot W function
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
close all 

C1 = 3.57;
C2 = 1.82;
C3 = 1.24;
ct = 3; 
% A2 = C1/ct;
% A3 = C2/ct;
A2 = C1;
A3 = C2;

syms theta 
A4_solve= 0:0.1:2;
for ind= 1:length(A4_solve)
B1_solve = 2*int(exp(-2*A2*theta/ct)*cos(A3*theta/ct+A4_solve(ind))*cos(A3*theta/ct),theta,0,inf);%alpha
B2_solve = 2*int(exp(-2*A2*theta/ct)*cos(A3*theta/ct+A4_solve(ind))*sin(A3*theta/ct),theta,0,inf); %beta
C3_solve(ind) = A4_solve(ind)+atan2(B2_solve,B1_solve);
end 

C3_solve = double(C3_solve);
[val,ind] = min(abs(C3_solve-C3));
A4 = A4_solve(ind);
figure()
plot(A4_solve,C3_solve)
hold on 
yline(C3)
hold off
xlabel('A4')
ylabel('phase shift with A4')
legend('phase shift using A4','C3')

% plot Falk method solution
B1_s = 2*int(exp(-2*A2*theta/ct)*cos(A3*theta/ct+A4)*cos(A3*theta/ct),theta,0,inf);%alpha
B2_s = 2*int(exp(-2*A2*theta/ct)*cos(A3*theta/ct+A4)*sin(A3*theta/ct),theta,0,inf); %beta
B1=double(B1_s);
B2=double(B2_s);
A1 = cos(A4)/(sqrt(cos(C3))*(B1^2+B2^2)^(1/4));


dx=0.01;
x=-4:dx:4;
W = (A1/cos(A4).*exp(-A2*abs(x)/ct).*cos(A3*abs(x)/ct+A4));
int_W=  sum(W.^2)*dx;


CX = (1/cos(C3)).*exp(-C1*abs(x)).*cos(C2*abs(x)+C3);

[acf,lags] = xcorr(W);
acf=acf/(sum(W.^2)*dx);


figure()
plot(x,W,'LineWidth',1.2)
xlabel('$\theta$','Interpreter','latex')
ylabel('W')

figure()
plot(lags*dx,acf)
hold on 
plot(x,CX)
hold off 
legend('auto of W','CX')
xlabel('$\Delta x$','Interpreter','latex')
ylabel('autocorrelation')







