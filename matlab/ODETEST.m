clear; clc;
fontSize = 20;
V0 = 25;
Vtip=7;
% max g =6 --> 58.86 m/s^2 lol
% envergadura = 8.84 --< R= 
tf = pi/2*4.42/Vtip;
syms V(t) x(t) z(t) mu(t) T(t) g Surf Cl(t) Cd(t) alpha rho m
eqs = [m*g*cos(mu(t)) + 1/2*rho*Surf*V(t)^2*Cl*alpha == 0,
       diff(x(t), t)==V(t)
       diff(z(t), t)==-0];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars,T(t), g, Surf, Cl(t), Cd(t), rho , m, alpha, mu(t));
m = 1e3;
g = 9.81;
Surf = 100;
alpha = 0.3;
mu = @(t) -Vtip/4.42*t+pi;
Cl = @(t) 0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3;
F = @(t,Y) F(t,Y,T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha, mu(t));
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [V0; 0; 0];
yp0 = [0; V0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
[t1, sol1]= ode15s(F, [t0, tf], y0, opt);
figure()
plot(t1, sol1(:,1)); hold on
plot(t1, sol1(:,2));
xlabel('Time')
legend('VSol', 'x(t)','Location','best')

top =9;
for i=1:top
    [Vfit_coff, S] = polyfit(t1, sol1(:,1), i);
    Svec(i)=sum(abs(S.R(:, i)));
    normr(i) = S.normr;
end
figure(); yyaxis left; plot(1:top, Svec); hold on; yyaxis right; plot(1:top, normr);
legend('$\Sum$ S', '|R^2|','Location','best', 'interpreter','latex')
title('Square root for phase 1')

figure()
plot(t1, sol1(:,1)); hold on;
[Vfit_coff] = polyfit(t1, sol1(:,1), 2);
P = poly2sym(Vfit_coff, t);
fplot(P);
xlabel('Convergence')
xlim([0 1]);

figure()
[Vfit_coff] = polyfit(sol1(:,2), sol1(:,1), 7);
P = poly2sym(Vfit_coff, t);
Cl = 0.12*t+0.33;
L = 1/2*rho*Surf*P^2*Cl*alpha;
fplot(L)
xlim([0 1]);
ylabel('Lift [N]')
xlabel('Time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
fontSize = 20;
V0 = 0;
Vtip=7;
% max g =6 --> 58.86 m/s^2 lol
% envergadura = 8.84 --< R= 
tfb = pi/2*4.42/Vtip;
syms V(t) x(t) z(t) mu(t) T(t) g Surf Cl(t) Cd(t) alpha rho m
eqs = [m*g*cos(mu(t)) - 1/2*rho*Surf*V(t)^2*Cl*alpha == 0,
       diff(x(t), t)==V(t)
       diff(z(t), t)==-0];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars,T(t), g, Surf, Cl(t), Cd(t), rho , m, alpha, mu(t));
m = 1e3;
g = 9.81;
Surf = 100;
alpha = 0.3;
mu = @(t) -Vtip/4.42*t+pi/2;
Cl = @(t) 0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3;
F = @(t,Y) F(t,Y,T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha, mu(t));
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [0; 0; 0];
yp0 = [1; V0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
[t1, sol1]= ode15s(F, [t0, tfb], y0, opt);
figure()
plot(t1, sol1(:,1)); hold on
plot(t1, sol1(:,2));
xlabel('Time')
legend('VSol', 'x(t)','Location','best')

top =9;
for i=1:top
    [Vfit_coff, S] = polyfit(t1, sol1(:,1), i);
    Svec(i)=sum(abs(S.R(:, i)));
    normr(i) = S.normr;
end
figure(); yyaxis left; plot(1:top, Svec); hold on; yyaxis right; plot(1:top, normr);
legend('$\Sum$ S', '|R^2|','Location','best', 'interpreter','latex')
title('Square root for phase 1')

figure()
plot(t1, sol1(:,1)); hold on;
[Vfit_coff] = polyfit(t1, sol1(:,1), 2);
P = poly2sym(Vfit_coff, t);
fplot(P);
xlabel('Convergence')
xlim([0 1]);

figure()
[Vfit_coff] = polyfit(sol1(:,2), sol1(:,1), 7);
P = poly2sym(Vfit_coff, t);
Cl = 0.12*t+0.33;
L = 1/2*rho*Surf*P^2*Cl*alpha;
fplot(L)
xlim([0 1]);
ylabel('Lift [N]')
xlabel('Time')