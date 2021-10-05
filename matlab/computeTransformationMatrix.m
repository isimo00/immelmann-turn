function R = computeTransformationMatrix(delta1,delta2,delta3)
R1 = rotationMatrix(delta1,3);
R2 = rotationMatrix(delta2,2);
R3 = rotationMatrix(delta3,1);
R = R3*R2*R1;
end

function R = rotationMatrix(angle,direction)
% Es crea una matriu R general que conté les combinacions necessàries per les
% matrius de rotació dels tres eixos; x és 1, y és 2 i z és 3.
% Les línies 18-20 fan aquesta conversií introduïnt 1 i 0 en R.
% Segons la direcció (1, 2 o 3) que es doni a la funció, aquesta
% reconverteix la matriu R en la matriu de rotació particular de l'eix en
% qüestió.
R = [cos(angle) sin(angle) -sin(angle);
    -sin(angle) cos(angle) sin(angle);
    sin(angle) -sin(angle) cos(angle)];
R(:,direction) = 0;
R(direction,:) = 0;
R(direction,direction) = 1;
end