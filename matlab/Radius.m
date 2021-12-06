function Radius
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
R = @(gamma) (2*6000*0.5.^2/(2*6000*9.81*cos(gamma)+1.225*20*0.5.^2*1.1));
fplot(R, [0 pi])
ylabel('Radius [m]')
xlabel('\gamma angle [rad]')
set(gca,'XTick',0:pi/4:pi)
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})

V = @(gamma) sqrt(2*6000*9.81*cos(gamma)*1/(2*6000*1-1.225*20*1.1));
fplot(V, [0 pi])
ylabel('V [m/s]')
xlabel('\gamma angle [rad]')
set(gca,'XTick',0:pi/4:pi)
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})

end