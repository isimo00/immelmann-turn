function [e11, e22] = intrinseques(e1,e2, mu)
%   Conversió d'equacions d'eixos wind a eixos horitzó
%   AUTHOR: Irene Simo Munoz 26/10/2021

e22 = simplify(expand(e1*sin(mu)+e2*cos(mu)));
e11 = simplify(expand(e1*cos(mu)-e2*sin(mu)));
end

