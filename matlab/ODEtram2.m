clear; clc;
fontSize = 20;
V0 = 16.5892;
R=6;
tf = R/V0 *pi;
syms V(t) gamma(t) x(t) z(t) T(t) g Surf Cl(t) Cd(t) alpha rho m
eqs = [T(t) - 1/2*rho*Surf*V(t)^2*Cd(t)*alpha - m*g*sin(gamma(t)) - 1e3*diff(V(t), t)==0,

       diff(x(t), t)==V(t)*cos(gamma(t)),
       diff(z(t), t)==-V(t)*sin(gamma(t))];
vars = [V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars, 'sparse', true);
F = odeFunction(F,vars, gamma(t), T(t), g, Surf, Cl(t), Cd(t), rho , m, alpha);
gamma = @(t) V0/R*t+(pi-V0/R*tf);
m = 1e3;
g = 9.81;
Surf = 50;
alpha = 0.3;
Cl = @(t) 0.12*alpha+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3+cos(t)*0.01;
F = @(t,Y) F(t,Y, gamma(t), T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [V0; 18.2108329349584; 0];
yp0 = [0; V0*cos(gamma(t0)); -V0*sin(gamma(t0))];
opt = odeset('mass',M,'InitialSlope',yp0);
[t2, sol2]= ode15s(F, [t0, tf], y0, opt);
figure()
plot(t2, sol2(:,1)); hold on
plot(t2, sol2(:,2)); hold on
plot(t2, sol2(:,3));
xlabel('Time')
legend('VSol','x(t)','-z(t)','Location','best')

top =9;
for i=1:top
    [Vfit_coff, S] = polyfit(t2, sol2(:,1), i);
    Svec(i)=sum(abs(S.R(:, i)));
    normr(i) = S.normr;
end
figure(); yyaxis left; plot(1:top, Svec); hold on; yyaxis right; plot(1:top, normr);

[Vfit_coff, S] = polyfit(t2, sol2(:,1), 5);
figure()
plot(t2, sol2(:,1)); hold on
P = poly2sym(Vfit_coff, t);
fplot(P);
xlim([0 tf]);


figure()
Cl = 0.12*t+0.33;
L = 1/2*rho*Surf*P^2*Cl*alpha;
fplot(L)
xlim([0 tf]);
ylabel('Lift [N]')
xlabel('Time')

figure()
Cd = 0.01+0.1*Cl^2;
D = 1/2*rho*Surf*P^2*Cd*alpha;
fplot(D)
xlim([0 tf]);
ylabel('Drag [N]')
xlabel('Time')

figure()
plot(sol2(:, 2), sol2(:, 3))
ylabel('-z(t)')
xlabel('x(t)')