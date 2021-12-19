V0 =70;
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
m = 1.3e3;
g = 9.81;
Surf = 13.69;
alpha = 0.3;
Cl = @(t) 0.12*alpha+0.33;
Cd = @(t) 0.01+0.1*Cl(t)^2;
rho = 1.225;
T = @(t) 1.5e3+cos(t)*0.01;
F = @(t,Y) F(t,Y, gamma(t), T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha);
t0 = 0;
y0 = [V0; 18.2108329349584; 0];
yp0 = [0; V0*cos(gamma(t0)); -V0*sin(gamma(t0))];
opt = odeset('mass',M,'InitialSlope',yp0);
[t2, sol2]= ode15s(F, [t0, tf], y0, opt);
figure()
plot(t2, sol2(:,1)); hold on
plot(t2, sol2(:,2)); hold on
plot(t2, abs(sol2(:,3)));xlabel('Time [s]'); ylabel('-')
legend('VSol','X(t)','-Z(t)','Location','best')

top =9;
for i=1:top
    [Vfit_coff, S] = polyfit(t2, sol2(:,1), i);
    Svec(i)=sum(abs(S.R(:, i)));
    normr(i) = S.normr;
end
figure(); yyaxis left; plot(1:top, Svec); ylabel('-'); hold on; 
yyaxis right; plot(1:top, normr); ylabel('-')
xlabel('Iterations');
legend('Sum of S','Squared regression','Location','best')

[Vfit_coff, S] = polyfit(t2, sol2(:,1), 5);
figure()
plot(t2, sol2(:,1), '-o'); hold on
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
plot(sol2(:, 2), sol2(:, 3))
ylabel('-Z(t)')
xlabel('X(t)')

%% Initial velocity effect
figure()
subplot(1, 2, 1)
for V0=70:0.1:70.4
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
m = 1.3e3;
g = 9.81;
Surf = 13.69;
alpha = 0.3;
Cl = @(t) 0.12*alpha+0.33;
Cd = @(t) 0.01+0.1*Cl(t)^2;
rho = 1.225;
T = @(t) 1.5e3+cos(t)*0.01;
F = @(t,Y) F(t,Y, gamma(t), T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha);
t0 = 0;
y0 = [V0; 18.2108329349584; 0];
yp0 = [0; V0*cos(gamma(t0)); -V0*sin(gamma(t0))];
opt = odeset('mass',M,'InitialSlope',yp0);
[t2, sol2]= ode15s(F, [t0, tf], y0, opt);
plot(t2, sol2(:,1)); hold on
end
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('V_0 =20 m/s','V_0 =20.1 m/s', 'V_0 =20.2 m/s', 'V_0 = 20.3 m/s', 'V_0 =20.4 m/s','Location','best')

subplot(1, 2, 2)
for R=6:1:10
V0=70;
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
m = 1.3e3;
g = 9.81;
Surf = 13.69;
alpha = 0.3;
Cl = @(t) 0.12*alpha+0.33;
Cd = @(t) 0.01+0.1*Cl(t)^2;
rho = 1.225;
T = @(t) 1.5e3+cos(t)*0.01;
F = @(t,Y) F(t,Y, gamma(t), T(t),m, g, Surf, Cl(t), Cd(t), rho, alpha);
t0 = 0;
y0 = [V0; 18.2108329349584; 0];
yp0 = [0; V0*cos(gamma(t0)); -V0*sin(gamma(t0))];
opt = odeset('mass',M,'InitialSlope',yp0);
[t2, sol2]= ode15s(F, [t0, tf], y0, opt);
plot(t2, sol2(:,1)); hold on
end
xlabel('Time')
legend('R =6 m','R =7 m', 'R =8 m', 'R =9 m', 'R =10 m','Location','best')