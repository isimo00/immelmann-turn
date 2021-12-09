function plot_trajectory
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% x =(@gamma) x;
% z =(@gamma) z;
% fplot(x, z);

x  = @(gamma) cos(gamma);
z =@(gamma) -sin(gamma);
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
hold off

x  = @(t) cos(t);
figure()
fplot(x, [pi/2 3*pi/2], 'LineWidth',1.5)
% set(gca,'YTick',pi/2:pi/4:3*pi/4)
% set(gca,'XTickLabel',{'\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2'})
ylabel('X horizon')
xlabel('\gamma')
hold off
z =@(gamma) -sin(gamma);
figure()
fplot(z, [pi/2 3*pi/2], 'LineWidth',1.5)
% set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2'})
ylabel('Z horizon')
xlabel('\gamma')
end