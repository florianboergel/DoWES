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