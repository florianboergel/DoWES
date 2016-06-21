import matplotlib
import pylab as pl
import numpy as np

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


for j in range(len(cp_list)):
    for i in range(len(cp_list[1][:])):
        writeData_cp.write(cp_list[j][i])
        writeData_thrust.write(thrust_list[j][i])
        writeData_torque.write(torque_list[j][i])

writeData_cp.close()
writeData_torque.close()
writeData_thrust.close()

cp = np.loadtxt(open("Output/cp.txt"))
thrust = np.loadtxt(open("Output/thrust.txt"))
torque = np.loadtxt(open("Output/torque.txt"))

pl.figure()
pl.plot(cp[0:38,0],cp[0:38,1])
pl.plot(cp[39:77,0],cp[39:77,1])
pl.plot(cp[78:116,0],cp[78:116,1])
pl.plot(cp[117:155,0],cp[117:155,1])
pl.plot(cp[156:195,0],cp[156:195,1])
pl.plot(cp[196:233,0],cp[196:233,1])
pl.plot(cp[234:272,0],cp[234:272,1])
axes = pl.gca()
axes.set_ylim([0,1])
pl.title('cp')
pl.savefig('cp.png')
pl.figure()
pl.plot(thrust[0:38,0],thrust[0:38,1])
pl.plot(thrust[39:77,0],thrust[39:77,1])
pl.plot(thrust[78:116,0],thrust[78:116,1])
pl.plot(thrust[117:155,0],thrust[117:155,1])
pl.plot(thrust[156:195,0],thrust[156:195,1])
pl.plot(thrust[196:233,0],thrust[196:233,1])
pl.plot(thrust[234:272,0],thrust[234:272,1])
axes = pl.gca()
pl.title('thrust')
pl.savefig('thrust.png')
pl.figure()
pl.plot(torque[0:38,0],torque[0:38,1])
pl.plot(torque[39:77,0],torque[39:77,1])
pl.plot(torque[78:116,0],torque[78:116,1])
pl.plot(torque[117:155,0],torque[117:155,1])
pl.plot(torque[156:195,0],torque[156:195,1])
pl.plot(torque[196:233,0],torque[196:233,1])
pl.plot(torque[234:272,0],torque[234:272,1])
axes = pl.gca()
pl.title('torque')
pl.savefig('torque.png')
pl.show()
