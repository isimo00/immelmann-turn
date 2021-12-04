% Irene Simo Munoz
% Immelmann turn
clear all; close; clc
paramEvolution
%% Symbolic definitions
T        = sym('Thrust', 'real');
W        = sym('Weight', 'real'); % Escalar, mòdul
Cl       = sym('Cl', 'real');
epsilon  = sym('epsilon', 'real');
nu       = sym('nu', 'real');
m        = sym('massa', 'real');
D        = sym('Drag', 'real');
Q        = sym('Q', 'real');
L        = sym('Lift', 'real');
V        = sym('Velocity', 'real');
Vdot     = sym('Vdot', 'real');
xi       = sym('xi', 'real');
xidot    = sym('xidot', 'real');
gamma    = sym('gamma', 'real');
gammadot = sym('gammadot', 'real');
mu       = sym('mu', 'real');
g        = sym('g', 'real');

rho      = sym('density', 'real');
S        = sym('Surface', 'real');  
Cd       = sym('Cd', 'real');
Cd_0     = sym('Cd_0', 'real');
k        = sym('k', 'real');

%% Maneuver study
% Equacions cinemàtiques i dinàmiques
% Condicions:
cond_vector1 = [gammadot, 0; % rectinineo
               gamma, -gamma% descenso
               Q, 0;        % vuelo simétrico
               T, 0;        % planeador
               xi, 0;       % plano vertical
               xidot, 0];   % plano vertical

% Equations
eq1 = T*cos(epsilon)*cos(nu)-D-m*g*sin(gamma)-m*Vdot; % ==0
eq2 = T*cos(epsilon)*sin(nu)-Q+m*g*cos(gamma)*sin(mu)+m*V*(gammadot*sin(mu)-xidot*cos(gamma)*cos(mu));
eq3 = -T*sin(epsilon)-L+m*g*cos(gamma)*cos(mu) + m*V*(gammadot*cos(mu)+xidot*cos(gamma)*sin(mu));

xe_dot = sym('xe_dot','real');
ye_dot = sym('ye_dot','real');
ze_dot = sym('ze_dot','real');

eq4 = xe_dot == V*cos(gamma)*cos(xi);
eq5 = ye_dot == V*cos(gamma)*sin(xi);
eq6 = ze_dot == -V*sin(gamma);

% eq_lift = L == 1/2*rho*V^2*S*Cl;
% eq_drag = D == 1/2*rho*V^2*S*Cd;
% eq_polar = Cd == Cd_0+kCl^2;
% eq_drag = subs(eq_drag, Cd, eq_polar);

eqv1 = [eq1; eq2; eq3; eq4; eq5; eq6];
% Substitució condicions en les equacions
for i=1:length(eqv1)
    for j=1:length(cond_vector1)
        eqv1(i)=subs(eqv1(i), [cond_vector1(j, 1) m*g m], [cond_vector1(j, 2) W W/g]);
    end
end
[eqv1(2), eqv1(3)]=wind2horizon(eqv1(2), eqv1(3), mu); % eixos vent a horitzó
%   Imposicio natural de mu=0 per planejador en vol rectilini (o pla
%   vertical?)
musol = solve(eqv1(2)==0, mu);
eqv1(2) = subs(eqv1(2), mu, musol);
eqv1(3) = subs(eqv1(3), mu, musol);
clearvars musol
fprintf('APARTAT 1:\n');
disp(eqv1);

dV = sym('dV', 'real');
dt = sym('dt', 'real');
eq = solve(eqv1(1), Vdot);
eq = eq == dV/dt;
eq = subs(eq, [D], [1/2*rho*V^2*S*(Cd_0+k*Cl^2)]);
LSol = solve(eqv1(3), L); % Warning solucions
eq = subs(eq, [Cl], [LSol/(0.5*rho*V^2*S)]);
int_temps1 = solve(eq, dt);
disp('dt=');
pretty(int_temps1)
fprintf('Entre V_basica i V_2.\n');
% temps1 = int(a, V); temps1 = subs(temps1, dV, 1);
% temps1 = simplify(expand(temps1));
% pretty(temps1)
disp('~~~~~~~~~~~~~~~~~~~~~~~~');

dh = sym('dh', 'real');
eq = lhs(eq) == -dV/dh*rhs(eqv1(6));
eq = subs(eq, [dh/dt], [-rhs(eqv1(6))]);
int_altura = solve(eq, dh);
disp('dh=')
pretty(int_altura);
fprintf('Entre V_basica i V_2.\n');
clearvars eq

fprintf('-----------------------\n-----------------------\n-----------------------\n');
