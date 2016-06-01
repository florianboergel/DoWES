import math as m
### constants ###
a_c = 0.2
# 1. Initialize a and a_0 to 0
a, a_ = 0
# 2. Evaluate the inflow angle alpha_a
alpha = m.atan((V_1*(1-a))/(rho*r*(1+a_)))
# 3. angle of attack alpha_a
alpha_a = alpha - (alpha_twist + alpha_pitch)
# 4. Interpolate the C_l,d - alpha_a curves at the angle of attack alpha_a
# 5. Compute thrust and torque coefficient C_thrust and C_torque
c_thrust = C_l*m.cos(alpha)+C_d*m.sin(alpha)
c_torque = C_l*m.sin(alpha)-C_d*m.cos(alpha)
# 6. Evaluate a and a_
sigma = N*c/(2*m.pi*r)
a_ = ((4*m.sin(alpha)*m.cos(alpha))/(sigma*c_torque)-1)
a_ = m.pow(a_,-1)
if a > a_c:
    a = 1/2*(2+K*(1-2*a_c)-m.sqrt(m.pow((K*(1-2*a_c)+2),2)+4*(K*m.pow(a_c,2)-1))
