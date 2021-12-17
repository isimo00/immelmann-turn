clear all; clc;
% V0 = 70;
% syms alpha(t) V(t) x(t) T(t) g S Cl(t) Cd(t) rho m
% eqs = [T(t)- 1/(2)*rho*S*V(t)^2*Cd*alpha(t) == diff(V(t),t),
%        m*g - 1/2*rho*S*V(t)^2*Cl*alpha(t) == 0,
%        diff(x(t), t)==V(t)];
% vars = [alpha(t) V(t) x(t)];
% [eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
% [M,F] = massMatrixForm(eqs,vars);
% M = odeFunction(M,vars, 'sparse', true);
% F = odeFunction(F,vars,T(t), g, S, Cl(t), Cd(t), rho , m)
% m = 1e3;
% g = 9.81;
% S = 100;
% Cl = @(t) 0.12*t+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
% Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
% rho = 1.225;
% T = @(t) 1.5e3+t*0.01;
% F = @(t,Y) F(t,Y,T(t),m, g, S, Cl(t), Cd(t), rho);
% t0 = 0;
% % y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
% y0 = [0; V0; 0];
% yp0 = [0; 0; 0];
% opt = odeset('mass',M,'InitialSlope',yp0);
% ode15s(F, [t0, 1], y0, opt)
% xlabel('Time')
% legend('VSol','\alpha', 'x(t)','Location','best')

syms x(t) y(t) z(t) T(t) b(t) c(t) S rho alpha m g
eqs = [T(t)- 1/2*rho*S*c(t)*x(t)^2*alpha - m*g*sin(y(t)) - diff(x(t), t)==0,
       cos(y(t)) - diff(y(t), t)== b(t)*x(t)^2,
       diff(x(t), t)==z(t)*cos(y(t))];
vars = [x(t) y(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars, T(t), b(t), c(t), rho, S, alpha, m, g);
T = @(t) 1.5e3+t*0.01;
b = @(t) 0.12*t+0.33;
c = @(t) 0.01+0.1*b(t)^2;
alpha =0.3;
rho=1.225;
S=100;
m = 1e3;
g=9.81;
F = @(t,Y) F(t,Y, T(t), b(t), c(t), rho, S, alpha, m, g);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [1; 0; 0];
yp0 = [0; 1; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
ode15s(F, [t0, 0.5], y0, opt);