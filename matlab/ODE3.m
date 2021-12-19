clear; clc;
fontSize = 20;
V0 = 70.194;
Vtip=7;
tf = pi*4.42/Vtip;
syms V(t) x(t) z(t) mu(t) T(t) g Surf Cl(t) Cd(t) alpha rho m
eqs = [m*g*cos(mu(t)) + 1/2*rho*Surf*V(t)^2*Cl*alpha == 0,
       diff(x(t), t)==V(t)
       diff(z(t), t)==-0];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars,T(t), g, Surf, Cl(t), Cd(t), rho , m, alpha, mu(t));
m = 1.3e3;
g = 9.81;
Surf = 13.69;
alpha = 0.3;
mu = @(t) -Vtip/4.42*t+pi;
Cl = @(t) 0.33;
Cd = @(t) 0.01+0.1*Cl(t)^2;
rho = 1.225;
T = @(t) 1.5e3;
F = @(t,Y) F(t,Y,T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha, mu(t));
t0 = 0;
y0 = [V0; 0; 0];
yp0 = [0; V0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
[t3, sol3]= ode15s(F, [t0, 1], y0, opt);
figure()
plot(t3, sol3(:,1)); hold on
plot(t3, sol3(:,2));
plot(t3, abs(sol3(:,3)));xlabel('Time [s]'); ylabel('-')
legend('VSol','X(t)','-Z(t)','Location','best')

top =9;
for i=1:top
    [Vfit_coff, S] = polyfit(t3, sol3(:,1), i);
    Svec(i)=sum(abs(S.R(:, i)));
    normr(i) = S.normr;
end
figure(); yyaxis left; plot(1:top, Svec); ylabel('-'); hold on; 
yyaxis right; plot(1:top, normr); ylabel('-')
xlabel('Iterations');
legend('Sum of S','Squared regression','Location','best')

figure()
plot(t3, sol3(:,1), '-o'); hold on;
[Vfit_coff] = polyfit(t3, sol3(:,1), 2);
P = poly2sym(Vfit_coff, t);
fplot(P, 'LineWidth', 2);
xlim([0 tf]);xlabel('Time [s]'); ylabel('Velocity [m/s]'); title('Curve fitting for the half loop');
legend('Set of ODE solved data','Fitted curve','Location','best')

figure()
subplot(1, 2, 1)
Cl = 0.12*t+0.33;
L = 1/2*rho*Surf*P^2*Cl*alpha;
fplot(L, 'LineWidth', 2)
xlim([0 tf]);
ylabel('Lift [N]')
xlabel('Time [s]')

subplot(1, 2, 2)
Cd = 0.01+0.1*Cl^2;
D = 1/2*rho*Surf*P^2*Cd*alpha;
fplot(D, 'LineWidth', 2)
xlim([0 tf]);
ylabel('Drag [N]')
xlabel('Time [s]')


figure()
plot(sol3(:, 2), sol3(:, 3))
ylabel('-Z(t)')
xlabel('X(t)')

figure()
fplot(cos(mu(t)), 'LineWidth', 2); hold on; xlim([0 tf]);xlabel('Time [s]');
fplot(sqrt(cos(mu(t))), 'Linewidth', 2)
legend('$cos(\mu)$','$\sqrt(cos(\mu))$','Location','best', 'interpreter', 'Latex')
