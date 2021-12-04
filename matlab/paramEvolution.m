function paramEvolution
% Evolution of different parameters
%   Detailed explanation goes here
%mu_val = @(t) pi;
syms gammaEV(t)
syms muEV(t)
gammaEV(t) = piecewise(0<t<0.75, 4/3*pi*t, 0.75<t, pi);
muEV(t) = piecewise(0<t<0.75, 0, 0.75<t, 4*pi*(t-0.75));
figure()
fplot(gammaEV, 'LineWidth',1.5); hold on
fplot(muEV, 'LineWidth',1.5);
xlim([0 1])
set(gca,'xtick',[])
xlabel('Time')
ylabel('Angle [rad]')
set(gca,'YTick',0:pi/4:pi)
set(gca,'YTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
legend('\gamma', '\mu', 'Location','northwest')

end