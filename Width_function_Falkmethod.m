% Bingchen Liu July 3
% Updated version of solving for Width function using Falk method and
% include scaling of ct in the W expression 


clear
close all 

ct=[2.5,3.5,4.5] ;%1:3;
for ind = 1:length(ct)
    [A1(ind),A2(ind),A3(ind),A4(ind),residue(ind)] = coef_solver(ct(ind));
end 


C1 = 3.57;
C2 = 1.82;
C3 = 1.24;
%% plot and validation 
% ================== plot solved W and check normalization 
close all 

for ind = 1:length(ct)
dx=0.01;
x=-5:dx:5;
x_nond(:,ind) = x/ct(ind);
W(:,ind) = (A1(ind)/cos(A4(ind)).*exp(-A2(ind)*abs(x)/ct(ind)).*cos(A3(ind)*abs(x)/ct(ind)+A4(ind)));
int_W(ind)=  sum(W(:,ind).^2)*dx;


% ================= calculate autocorrelation and compare with CX 
x_lag = -5:dx:5;
x_lag_nond(:,ind) = x_lag/ct(ind);
CX(:,ind) = (1/cos(C3)).*exp(-C1*abs(x_lag)/ct(ind)).*cos(C2*abs(x_lag)/ct(ind)+C3);

[acf(:,ind),lags(:,ind)] = xcorr(W(:,ind));
acf(:,ind)=acf(:,ind)/(sum(W(:,ind).^2)); %doesn't *dx due to the xcorr normalization

end 



%% plot 


figure()
for ind = 1:length(ct)
plot(x,W(:,ind))
xlabel('x(m) -- dimensional')
ylabel('W')
hold on 
end 
hold off 



figure()
for ind = 1:2
subplot(211)
plot(lags(:,ind)*dx,acf(:,ind),'LineStyle','-','LineWidth',1.2)
hold on 
plot(x_lag,CX(:,ind),'LineStyle','--','LineWidth',1.2)
hold off 
xlabel('$\Delta x$','Interpreter','latex')
ylabel('autocorrelation')
hold on 
subplot(212)
plot(lags(:,ind)*dx/ct(ind),acf(:,ind),'LineStyle','-','LineWidth',1.2)
hold on 
plot(x_lag_nond(:,ind),CX(:,ind),'LineStyle','--','LineWidth',1.2)
hold off 
xlabel('$\Delta x/ct$','Interpreter','latex')
ylabel('autocorrelation')
hold on 
end 
hold off 


%% function 

function [A1,A2,A3,A4,residue] = coef_solver(ct)

C1 = 3.57;
C2 = 1.82;
C3 = 1.24;
% A2 = C1/ct;
% A3 = C2/ct;
A2 = C1;
A3 = C2;

% =============== solve for A4 numerically by calculating A4+arctan(B2/B1)
% and compare with C3
syms theta 
A4_solve= 0:0.1:2;
for ind= 1:length(A4_solve)
B1_solve = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4_solve(ind))*cos(A3*theta),theta,0,inf); %note this is nondimensional theta
B2_solve = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4_solve(ind))*sin(A3*theta),theta,0,inf); 
C3_solve(ind) = A4_solve(ind)+atan2(B2_solve,B1_solve);
end 

C3_solve = double(C3_solve);
[residue,ind] = min(abs(C3_solve-C3));
A4 = A4_solve(ind);


% ============== solve for A1 using normalization and calculate in B1 and B2
% using A4

B1_s = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4)*cos(A3*theta),theta,0,inf);%alpha
B2_s = 2*int(exp(-2*A2*theta)*cos(A3*theta+A4)*sin(A3*theta),theta,0,inf); %beta
B1=double(B1_s);
B2=double(B2_s);
A1 = cos(A4)/(sqrt(cos(C3))*(B1^2+B2^2)^(1/4)*sqrt(ct));
end 

