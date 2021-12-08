function trejectory(x,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% x =(@gamma) x;
% z =(@gamma) z;
% fplot(x, z);

x  = @(gamma) cos(gamma);
z =@(gamma) sin(gamma);
figure()
fplot(x, z, [pi/2 3*pi/2], 'LineWidth',1.5);
hold on
x  = @(gamma) 1*gamma;
z =@(gamma) 1;
fplot(x, z, [0 pi/4], 'LineWidth',1.5);
hold on
x  = @(gamma) 1*gamma;
z =@(gamma) -1;
fplot(x, z, [0 pi/4], 'LineWidth',1.5);
xlim([-1.1 1.1])
ylim([-1.1 1.1])
set(gca,'ytick',[])
set(gca,'xtick',[])
xlabel('X horizon')
ylabel('-Z horizon')
end