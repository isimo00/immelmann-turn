function paramEvolution
% Author: Irene Simo Munoz
% Date: Fall 2021
% Evolution of different parameters trhoughtout time
%   This function requires no input and produces no output to the global
%   code. It creates a figure with the resulting evolution.

syms gammaEV(t)
syms muEV(t)
gammaEV(t) = piecewise(0<t<0.75, 4/3*pi*t, 0.75<t, pi);
muEV(t) = piecewise(0<t<0.75, pi, 0.75<t, -4*pi*(t-0.75)+pi);
figure()
fplot(gammaEV, 'LineWidth',1.5); hold on
fplot(muEV, 'LineWidth',1.5);
xlim([0 1])
set(gca,'xtick',[])
xlabel('Time')
ylabel('Angle [rad]')
set(gca,'YTick',0:pi/4:pi)
set(gca,'YTickLabel',{'0','\pi/4','\pi/2','3\pi/2','\pi'})
legend('\gamma', '\mu', 'Location','northwest')
hold off

end