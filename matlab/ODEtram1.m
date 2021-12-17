clear all; clc;
V0 = 70;
syms V(t) x(t) z(t) T(t) g S Cl(t) Cd(t) alpha rho m
eqs = [m*g - 1/2*rho*S*V(t)^2*Cl*alpha == 0,
       diff(x(t), t)==V(t)
       diff(z(t), t)==-0];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars,T(t), g, S, Cl(t), Cd(t), rho , m, alpha);
m = 1e3;
g = 9.81;
S = 100;
alpha = 0.3;
Cl = @(t) 0.12*t+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3+t*0.01;
F = @(t,Y) F(t,Y,T(t),m, g, S, Cl(t), Cd(t), rho, alpha);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [V0; 0; 0];
yp0 = [0; V0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
[t1, sol1]= ode15s(F, [t0, 1], y0, opt);
figure()
plot(t1, sol1(:,1)); hold on
plot(t1, sol1(:,2));
xlabel('Time')
legend('VSol', 'x(t)','Location','best')

figure()
[Vfit_coff] = polyfit(sol1(:,2), sol1(:,1), 7);
P = poly2sym(Vfit_coff, t);
Cl = 0.12*t+0.33;
L = 1/2*rho*S*P^2*Cl*alpha;
fplot(L)
xlim([0 1]);
ylabel('Lift [N]')
xlabel('Time')

figure()
Cd = 0.01+0.1*Cl^2;
D = 1/2*rho*S*P^2*Cd*alpha;
fplot(D)
xlim([0 1]);
ylabel('Drag [N]')
xlabel('Time')

figure()
plot(sol1(:, 2), sol1(:, 3))
ylabel('-z(t)')
xlabel('x(t)')