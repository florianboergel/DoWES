a = 0;
a' = 0;
rotorRadius = 50;
numBlades = 3;
v1 = 10;
airDensity = 1.225;
alpha_twist = 0;
alpha_pitch = 0;
alpha = atan((v1(1-a))/(omega*r*(1+a')));

angleOfAttack = alpha - (alpha_twist + alpha_pitch)