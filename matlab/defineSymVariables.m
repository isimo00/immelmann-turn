function [outputArg1,outputArg2] = defineSymVariables
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Definicions per variables de les equacions
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
end