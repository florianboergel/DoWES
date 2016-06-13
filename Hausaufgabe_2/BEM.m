rotor_diameter = 100
number_of_blade = 3
ele_conv_eff = 0.95
design_ws = 10
design_tip_speed_ratio = 7
design_lift_coefficient = 0.65
air_density = 1.225
a_0 = 0
inflow_angle = 0
angle_of_attack = 0

inflow_angle = atan(design_ws*(1-a_0)/(design_tip_speed_ratio*0.5*rotor_diameter))
