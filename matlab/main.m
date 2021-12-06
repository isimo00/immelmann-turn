% Author: Irene Simo Munoz
% Date: Fall 2021
% Project: Immelmann turn flight mechanics study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close; clc
%% 
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
% Cinematic and dynamic equations
% Equations
eq1 = T*cos(epsilon)*cos(nu)-D-W*sin(gamma)-W/g*Vdot; % ==0
eq2 = T*cos(epsilon)*sin(nu)-Q+W*cos(gamma)*sin(mu)+W/g*V*(gammadot*sin(mu)-xidot*cos(gamma)*cos(mu));
eq3 = -T*sin(epsilon)-L+W*cos(gamma)*cos(mu) + W/g*V*(gammadot*cos(mu)+xidot*cos(gamma)*sin(mu));

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

% Conditions:

eqv = [eq1; eq2; eq3; eq4; eq5; eq6];
% Substitució condicions en les equacions
%% First phase- cruise
cond_vector1 = [Q, 0;       % vuelo simétrico
               nu, 0;       % vuelo simétrico
               gamma, 0     % cruise
               gammadot, 0; % cruise
               mu, 0;       % levelled flight
               epsilon, 0;  % thrust parallel to x wind
               xi, 0;       % plano vertical
               xidot, 0];   % plano vertical
% cond_vector1a = [Q, 0;       % vuelo simétrico
%                nu, 0;       % vuelo simétrico
%                gamma, 0        % initial straig flight
%                mu, 0;
%                epsilon, 0;  % thrust parallel to x wind
%                xi, 0;       % plano vertical
%                xidot, 0];   % plano vertical
% cond_vector1b = [Q, 0;       % vuelo simétrico
%                nu, 0;       % vuelo simétrico
%                mu, 0;
%                gamma, pi/2;
%                epsilon, 0;   % thrust parallel to x wind
%                xi, 0;       % plano vertical
%                xidot, 0];   % plano vertical
% 
% cond_vector1c = [Q, 0;       % vuelo simétrico
%                nu, 0;       % vuelo simétrico
%                mu, pi;       % inverted flight
%                gamma, pi;
%                epsilon, 0;  % thrust parallel to x wind
%                xi, 0;       % plano vertical
%                xidot, 0];   % plano vertical
[eqv1] = studyTheseConditions(eqv,cond_vector1, mu);
fprintf('APARTAT 1:\n');
disp(eqv1);

% Evolution

%% Second phase - semicircle
cond_vector2 = [Q, 0;       % vuelo simétrico
               nu, 0;       % vuelo simétrico
               mu, 0;       % levelled flight
               epsilon, 0;  % thrust parallel to x wind
               xi, 0;       % plano vertical
               xidot, 0];   % plano vertical
[eqv2] = studyTheseConditions(eqv,cond_vector2, mu);
fprintf('APARTAT 2:\n');
disp(eqv2);

% lift_centripletal= isolate(eqv2(3)==0, L)
% syms R
% syms Cl
% R = @(gamma) 2*m*V^2/(2*m*g*cos(gamma)+rho*S*V^2*Cl);



%% Third phase - turnaround
cond_vector3 = [Q, 0;       % vuelo simétrico
               nu, 0;       % vuelo simétrico
               mu, -pi;       % levelled flight
               gamma, 0;
               gammadot, 0;
               epsilon, 0;  % thrust parallel to x wind
               xi, 0;       % plano vertical
               xidot, 0];   % plano vertical

[eqv3] = studyTheseConditions(eqv,cond_vector3, mu);
fprintf('APARTAT 3:\n');
disp(eqv3)
%% Forth phase - cruise
cond_vector4 = [Q, 0;       % vuelo simétrico
               nu, 0;       % vuelo simétrico
               mu, 0;       % levelled flight
               epsilon, pi;  % thrust parallel to x wind
               gamma, 0;
               gammadot, 0;
               xi, 0;       % plano vertical
               xidot, 0];   % plano vertical
[eqv4] = studyTheseConditions(eqv,cond_vector4, mu);
fprintf('APARTAT 4:\n');
disp(eqv4)
%%

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
