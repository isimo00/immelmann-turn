function [eqv] = studyTheseConditions(eqv,cond_vector, mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(eqv)
    for j=1:length(cond_vector)
        eqv(i)=subs(eqv(i), [cond_vector(j, 1)], [cond_vector(j, 2)]);
    end
end
%[eqv(2), eqv(3)]=wind2horizon(eqv(2), eqv(3), mu); % eixos vent a horitzó
end