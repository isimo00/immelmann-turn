function Radius
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure()
R = @(V) (2*1.3e3*V^2/(2*1.3e3*9.81+1.225*13.69*V^2*0.12));
subplot(1,2,1);
fplot(R, 'LineWidth',2'); xlim([0 25])
ylabel('Radius [m]')
xlabel('Velocity [m/s]')
set(gca,'YTickLabel',{})
set(gca,'XTickLabel',{})
set(gca,'XTick',0)
set(gca,'YTick',0)

R = @(alpha) (2*1.3e3*10^2/(2*1.3e3*9.81+1.225*13.69*10^2*(0.12+0.33*alpha)));
subplot(1,2,2);
fplot(R, [0 pi], 'LineWidth',2');
xlabel('\alpha angle [rad]')
set(gca,'XTick',0)
set(gca,'YTick',0)
set(gca,'YTickLabel',{})
set(gca,'XTickLabel',{})

end