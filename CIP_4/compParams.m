%CIP 4 Task c
set(0,'DefaultTextInterpreter', 'latex');

E=211000000000;
D=5;
l=100;
mTop=323000;
rho=7850;
t=1;
rrs=12.59/60;
f0=rrs*1.1;

func = @(t) sqrt(3*E*pi*D^3*t/(l^3*8*(mTop+0.25*rho*pi*D+t*l)))-0.23081*2*pi;
t = fsolve(func,t);

mTow= rho*pi*D*t*l;
towerCost = mTow/1000*500;

%plotting f_0 as function of t
figure;
hold on;
fplot( @(t) sqrt(3*E*pi*D^3*t/(l^3*8*(mTop+0.25*rho*pi*D+t*l)))/(2*pi), [0 0.5]);
plot(t,f0,'*');
xlabel('Wall thickness t [m]','FontSize',12);
ylabel('Eigenfrequency f0 [Hz]','FontSize',12);
title('Relationship: Eigenfrequency-Thickness', 'FontSize',14);
saveas(gcf,strcat('figures/eigenfrequency.jpg'));
hold off;

%campbell diagram
figure();
hold on;
fplot( @(x) f0, [0 0.5]);
fplot( @(x) 1*x, [0 0.5]);
fplot( @(x) 2*x, [0 0.5]);
fplot( @(x) 3*x, [0 0.5]);
legend('Tower bending','\Omega','2\Omega','3\Omega','location','Northwest');
xlabel('Rotor Speed [Hz]','FontSize',12);
ylabel('Eigenfrequencies $\Omega$ [Hz]','FontSize',12);
title('Campbell Diagram', 'FontSize',14);
ax2 = axes('Position',[0 0.1 0.8 0.001],'Color','none')

saveas(gcf,strcat('figures/campbell.jpg'));
hold off;

