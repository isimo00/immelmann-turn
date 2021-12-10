clear all; clc;
% Ode solver test
syms gamm(t) V(t)
m= 3e4;
T=2e5;
D=1e5;
g=9.81;
L=2e5;

% ode1 = diff(V) == 1/m*(T - D -m*g*sin(gamm));
% ode2 = diff(gamm) == (L-m*g*cos(gamm))/(m*V);
% odes = [ode1; ode2];
% % % cond1 = u(0) == 0;
% % % cond2 = v(0) == 1;
% % % conds = [cond1; cond2];
% [uSol(t),vSol(t)] = dsolve(odes, conds);
% %[uSol(t),vSol(t)] = dsolve(odes);
% fplot(uSol)
% hold on
% fplot(vSol)
% grid on
% legend('uSol','vSol','Location','best')

% syms x1(t) x2(t) T D L(t)
% eqs = [diff(V(t),t) == 1/m*(T - D -m*g*sin(gamm)),...
%        diff(gamm(t), t) == (L(t)-m*g*cos(gamm))/(m*V)];
% vars = [V(t) gamm(t)];
% [M,F] = massMatrixForm(eqs,vars)
% M = odeFunction(M,vars);
% F = odeFunction(F,vars, T, D, L(t));
% T = 2e5;
% D = 1e5;
% L = @(t) cos(t)/(1+t^2);
% F = @(t,Y) F(t,Y,T,D,L(t));
% t0 = 0;
% y0 = [-L(t0)*sin(0.1); L(t0)*cos(0.1)];
% % yp0 = [T*y0(1) + D*y0(2)^2; 1.234];
% % opt = odeset('mass',M,'InitialSlope',yp0);
% % ode15s(F, [t0, 1], y0, opt)
% ode15s(F, [t0, 1], y0)
% legend('VSol','\gammaSol','Location','best')
% xlabel('Time');
% ylabel('?')

syms x1(t) x2(t) T rho S
eqs = [diff(V(t),t) == 1/m*(T - 0.5*rho*S*V(t)^2*0.005 -m*g*sin(gamm)),...
       diff(gamm(t), t) == (0.5*rho*S*V(t)^2*0.04-m*g*cos(gamm))/(m*V)];
vars = [V(t) gamm(t)];
[M,F] = massMatrixForm(eqs,vars)
M = odeFunction(M,vars);
F = odeFunction(F,vars, T, rho, S);
T = 2e5;
rho = 1.225;
S = 20;
F = @(t,Y) F(t,Y,T, rho, S);
t0 = 0;
% y0 = [-L(t0)*sin(0.1); L(t0)*cos(0.1)];
y0 = [sin(0.1); cos(0.1)];
% yp0 = [T*y0(1) + D*y0(2)^2; 1.234];
% opt = odeset('mass',M,'InitialSlope',yp0);
% ode15s(F, [t0, 1], y0, opt)
ode15s(F, [t0, 1], y0)
legend('VSol','\gammaSol','Location','best')
xlabel('Time');
ylabel('?')