clear all; clc;
V0 = 70;
syms alpha(t) V(t) x(t) T(t) g S Cl(t) Cd(t) rho m
eqs = [T(t)- 1/(2)*rho*S*V(t)^2*Cd*alpha(t) == diff(V(t),t),
       m*g - 1/2*rho*S*V(t)^2*Cl*alpha(t) == 0,
       diff(x(t), t)==V(t)];
vars = [alpha(t) V(t) x(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars);
F = odeFunction(F,vars,T(t), g, S, Cl(t), Cd(t), rho , m)
m = 1e3;
g = 9.81;
S = 100;
Cl = @(t) 0.12*t+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3+t*0.01;
F = @(t,Y) F(t,Y,T(t),m, g, S, Cl(t), Cd(t), rho);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [0; V0; 0];
yp0 = [0; 0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
ode15s(F, [t0, 1], y0, opt)
xlabel('Time')
legend('VSol','\alpha', 'x(t)','Location','best')

R = 3;
syms alpha(t) V(t) x(t) z(t) T(t) g S Cl(t) Cd(t) rho m
eqs = [T(t)- 1/(2)*rho*S*V(t)^2*Cd*alpha(t)-g*sin(V0/R) == diff(V(t),t),
       m*g - 1/2*rho*S*V(t)^2*Cl*alpha(t) +m*V(t)*0 == 0,
       diff(x(t), t)==cos(V0/R)*V(t),
       diff(z(t), t)==sin(V0/R)*V(t)];
vars = [alpha(t) V(t) x(t) z(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars);
F = odeFunction(F,vars,T(t), g, S, Cl(t), Cd(t), rho , m)
m = 1e3;
g = 9.81;
S = 100;
Cl = @(t) 0.12*t+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3+t*0.01;
F = @(t,Y) F(t,Y,T(t),m, g, S, Cl(t), Cd(t), rho);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [0; V0; 0];
yp0 = [0; 0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
ode15s(F, [t0, 1], y0, opt)
xlabel('Time')
legend('VSol','\alpha', 'x(t)','Location','best')
%%
eq1 = equacions;

% Variables que interessen
syms xe1(t) ze1(t) gamma_1(t) T1(t);

% Se substitueixen valors
eq1 = subs(eq1, [xi xi_dot], [0 0]); %no hi ha guinyada ni canvi en ella
eq1 = subs(eq1, [Q mu],[0 0]); %vol simetric, ales a nivell

% Es força una trajectoria circular amb radi R, V constant
eq1 = subs(eq1, [gamma_dot, V_dot], [V0/R, 0]);
eq1 = subs(eq1, V, V0);
V0 = 70;
syms alpha(t) V(t) x(t) T(t) g S Cl(t) Cd(t) rho m
eqs = [T(t)- 1/(2)*rho*S*V(t)^2*Cd*alpha(t) == diff(V(t),t),
       m*g - 1/2*rho*S*V(t)^2*Cl*alpha(t) == 0,
       diff(x(t), t)==V(t)];
vars = [alpha(t) V(t) x(t)];
[eqs, vars] = reduceDifferentialOrder(eqs, vars); % cAL?
[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars);
F = odeFunction(F,vars,T(t), g, S, Cl(t), Cd(t), rho , m)
m = 1e3;
g = 9.81;
S = 100;
Cl = @(t) 0.12*t+0.33; % https://en.wikipedia.org/wiki/ENAER_T-35_Pillán
Cd = @(t) 0.01+0.1*Cl(t)^2; % http://airfoiltools.com/airfoil/details?airfoil=naca652415-il
rho = 1.225;
T = @(t) 1.5e3+t*0.01;
F = @(t,Y) F(t,Y,T(t),m, g, S, Cl(t), Cd(t), rho);
t0 = 0;
% y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
y0 = [0; V0; 0];
yp0 = [0; 0; 0];
opt = odeset('mass',M,'InitialSlope',yp0);
ode15s(F, [t0, 1], y0, opt)
xlabel('Time')
legend('VSol','\alpha', 'x(t)','Location','best')
% Valors raonables
eq1(1) = valors_raonables(eq1(1), false);
eq1(3) = valors_raonables(eq1(3), false);

% Valor d'alpha
alpha1eq = solve(eq1(3), alpha);
eq1(1) = subs(eq1(1), alpha, alpha1eq);

% De sym a equcions diferencials
edos1(1) = subs(eq1(1), [gamma_ T], [gamma_1 T1]);
edos1(2) = diff(xe1) == subs(solve(eq1(4), xe_dot), gamma_, gamma_1);
edos1(3) = diff(ze1) == subs(solve(eq1(6), ze_dot), gamma_, gamma_1);
edos1(4) = xe1 == R*sin(gamma_1);

vars1 = [xe1 ze1 gamma_1 T1];
[edos1, vars1] = reduceDifferentialOrder(edos1, vars1);
[M1, F1] = massMatrixForm(edos1,vars1);
M1 = odeFunction(M1, vars1, 'sparse', true);
F1 = odeFunction(F1, vars1);

% Condicions de contorn i solucio del sistema
T0 = solve(eq1(1), T);
T0 = subs(T0, gamma_, 0);

diffT0 = subs(diff(edos1(1)), diff(gamma_1, t), V0/R);
diffT0 = subs(diffT0, gamma_1, 0);
diffT0 = subs(diffT0, diff(T1, t), T_dot);
diffT0 = solve(diffT0, T_dot);

y_0 = [x0; z0; 0; T0];
y_0 = double(y_0);

yp0 = zeros(size(vars1));
yp0(1) = V0;
yp0(3) = V0/R;
yp0(4) = diffT0;
yp0 = double(yp0);

tf1 = 0.9425; %el temps que es tarda en arribar a gamma=pi/2
tspan1 = [t0 tf1];
opt1 = odeset('mass',M1, 'InitialSlope',yp0);
[t1, sol1] = ode15s(F1, tspan1, y_0, opt1);

% Plots

% figure;
% plot(t1, sol1(:,1));
% grid;
% title('Tram 1')
% xlabel('Temps (s)');
% ylabel('Posició Xe (m)');

% figure;
% plot(t1, sol1(:,2));
% grid;
% title('Tram 1')
% xlabel('Temps (s)');
% ylabel('Posició Ze (m)');
% set(gca, 'YDir','reverse');

% figure;
% plot(sol1(:,1), sol1(:,2));
% grid;
% title('Tram 1');
% xlabel('Posició Xe (m)');
% ylabel('Posició Ze(m)');
% set(gca, 'YDir','reverse');

% figure;
% plot(t1, sol1(:,4));
% grid;
% title('Tram 1')
% xlabel('Temps (s)');
% ylabel('Thrust (N)');

L1 = zeros(size(sol1(:,1)));
D1 = zeros(size(sol1(:,1)));
alpha1 = zeros(size(sol1(:,1)));
for i=1:size(sol1(:,1))
    alpha1(i) = subs(alpha1eq, gamma_, sol1(i,3));
    L1(i) = subs(L_eq, [V, alpha], [V0, alpha1(i)]);
    D1(i) = subs(D_eq, [V, alpha], [V0, alpha1(i)]);
end

% figure;
% plot(t1, sol1(:,3));
% grid;
% hold on;
% plot(t1, alpha1);
% title('Tram 1')
% xlabel('Temps (s)');
% ylabel('angle (rad)');
% legend('Gamma', 'Alpha', 'location', 'best');

% figure;
% plot(t1, D1);
% title('Tram 1');
% grid;
% xlabel('Temps (s)');
% ylabel('Drag (N)');

% figure;
% plot(t1, L1);
% grid;
% title('Tram 1');
% xlabel('Temps (s)');
% ylabel('Lift (N)');

