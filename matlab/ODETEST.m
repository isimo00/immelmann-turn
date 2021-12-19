clear; clc;
V0 = 25;
syms V(t) x(t) z(t) T(t) Cl(t) Cd(t) m alpha rho g Surf
eqs = [m*g - 1/2*rho*Surf*V(t)^2*Cl*alpha == 0,
       diff(x(t), t)==V(t)
       diff(z(t), t)==-0];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars,T(t), g, Surf, Cl(t), Cd(t), rho , m, alpha);
m = 1.3e3;
g = 9.81;
Surf = 13.69;
alpha = 0.3;
Cl = @(t) 0.12*alpha+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3;
F = @(t,Y) F(t,Y,T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha);
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
Vfit = poly2sym(Vfit_coff, t);
fplot(Vfit);
ylabel('Velocity [m/s]')
xlabel('Time [s]')
title('Convergence')
xlim([0 1]);

figure()
Cl = 0.12;
L = 1/2*rho*Surf*Vfit^2*Cl*alpha;
fplot(L, 'LineWidth',2)
xlim([0 1]);
ylabel('Lift [N]')
xlabel('Time [s]')

figure()
Cd = 0.01+0.1*Cl^2;
D = 1/2*rho*Surf*Vfit^2*Cd*alpha;
fplot(D, 'LineWidth',2)
xlim([0 1]);
ylabel('Drag [N]')
xlabel('Time [s]')

figure()
plot(sol1(:, 2), sol1(:, 3))
ylabel('-z(t)')
xlabel('x(t)')