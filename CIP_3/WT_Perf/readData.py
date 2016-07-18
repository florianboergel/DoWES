import matplotlib
import pylab as pl
import numpy as np
from scipy.optimize import fsolve

rho = 1.225
A = np.pi * 54.00 * 54.00
Ud = 11

file = open('WT_test.oup')
data = file.readlines()

writeData_cp = open('Output/cp.txt', 'w')
writeData_thrust = open('Output/thrust.txt', 'w')
writeData_torque = open('Output/torque.txt', 'w')


cp_list = []
torque_list = []
thrust_list = []

# Lines depending on inputfile
for i in range(0,7):
    cp_list.append(data[20+i*145:59+i*145])
    thrust_list.append(data[112+i*145:151+i*145])
    torque_list.append(data[66+i*145:105+i*145])

# Extract data to output file
for j in range(len(cp_list)):
    for i in range(len(cp_list[1][:])):
        writeData_cp.write(cp_list[j][i])
        writeData_thrust.write(thrust_list[j][i])
        writeData_torque.write(torque_list[j][i])

writeData_cp.close()
writeData_torque.close()
writeData_thrust.close()

# calc and plot
cp = np.loadtxt(open("Output/cp.txt"))
thrust = np.loadtxt(open("Output/thrust.txt"))
torque = np.loadtxt(open("Output/torque.txt"))


for i in range(len(thrust[:,0])):
	#a_p = fsolve(func, 0.4)
	thrust[i,1] = 1000*2*thrust[i,1]*(thrust[i,0])**2/(rho*54.00**4*np.pi*(2*np.pi*14.50/60)**2)
	torque[i,1] = 1000*2*torque[i,1]*(torque[i,0])**2/(rho*54.00**5*np.pi*(2*np.pi*14.50/60)**2)


pl.figure()
pl.plot(cp[0:38,0],cp[0:38,1], label='0 deg.')
pl.plot(cp[39:77,0],cp[39:77,1], label='5 deg.')
pl.plot(cp[78:116,0],cp[78:116,1], label='10 deg.')
pl.plot(cp[117:155,0],cp[117:155,1], label='15 deg.')
pl.plot(cp[156:195,0],cp[156:195,1], label='20 deg.')
pl.plot(cp[234:272,0],cp[234:272,1], label='30 deg.')
pl.legend()
pl.xlabel('Tipspeed Ratio \lambda')
pl.ylabel('C_p')
axes = pl.gca()
axes.set_ylim([0,1])
pl.title('Power coefficient')
pl.savefig('Output/cp.png')
pl.figure()
pl.plot(thrust[0:38,0],thrust[0:38,1], label='0 deg.')
pl.plot(thrust[39:77,0],thrust[39:77,1], label='5 deg.')
pl.plot(thrust[78:116,0],thrust[78:116,1], label='10 deg.')
pl.plot(thrust[117:155,0],thrust[117:155,1], label='15 deg.')
pl.plot(thrust[156:195,0],thrust[156:195,1], label='20 deg.')
pl.plot(thrust[234:272,0],thrust[234:272,1], label='30 deg.')
axes = pl.gca()
axes.set_ylim([0,4])
pl.title('Thrust coefficient')
pl.legend()
pl.xlabel('Tipspeed Ratio \lambda')
pl.ylabel('C_t')
pl.savefig('Output/thrust.png')
pl.figure()
pl.plot(torque[0:38,0],torque[0:38,1], label='0 deg.')
pl.plot(torque[39:77,0],torque[39:77,1], label='5 deg.')
pl.plot(torque[78:116,0],torque[78:116,1], label='10 deg.')
pl.plot(torque[117:155,0],torque[117:155,1], label='15 deg.')
pl.plot(torque[156:195,0],torque[156:195,1], label='20 deg.')
pl.plot(torque[234:272,0],torque[234:272,1], label='30 deg.')
axes = pl.gca()
pl.title('Torque coefficient')
pl.legend()
pl.xlabel('Tipspeed Ratio \lambda')
pl.ylabel('C_q')
axes.set_ylim([0,0.08])
pl.savefig('Output/torque.png')

# Task 8
pl.figure()
pl.plot(cp[0:38,0],cp[0:38,1],label='0 deg.')
pl.plot(cp[39:77,0],cp[39:77,1],label='5 deg.')
pl.plot(cp[78:116,0],cp[78:116,1],label='10 deg.')
pl.plot(cp[117:155,0],cp[117:155,1],label='15 deg.')
pl.plot(cp[156:195,0],cp[156:195,1],label='20 deg.')
pl.plot(cp[196:233,0],cp[196:233,1],label='25 deg.')
pl.plot(cp[234:272,0],cp[234:272,1],label='30 deg.')
pl.axhline(y=0.272, xmin=0, xmax=1, hold=None, label='cp = 0.272')
pl.legend()
axes = pl.gca()
axes.set_ylim([0,1])
pl.title('Validate calc cp')
pl.savefig('Output/validated_cp.png')
#pl.show()


# Task 5

thrust = 392.262
torque = 968.677
RPM = 11.86
cp = 0.543895
tsp = 8.2
c_t = 0.0

2054.846

c_t = 2*thrust*tsp**2*1000/(rho*54.00**4*np.pi*(2*np.pi*RPM/60)**2)
c_q = 2*torque*tsp**2*1000/(rho*54.00**5*np.pi*(2*np.pi*RPM/60)**2)


theo_pow = 0.5*16/27*rho*np.pi*54.00**2*13**3
actual = 3500000/(0.5*rho*np.pi*54.00**2*13**3)
rpm_310 = 5 * 13 /(54.00*2*np.pi) * 60

print c_t, c_q, theo_pow, actual,rpm_310
