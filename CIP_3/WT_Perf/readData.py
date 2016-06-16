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
pl.plot(cp[39:39*2,0],cp[39:39*2,1])
pl.title('cp')
pl.figure()
pl.plot(thrust[0:38,0],thrust[0:38,1])
pl.plot(thrust[39:39*2,0],thrust[39:39*2,1])
pl.title('thrust')
pl.figure()
pl.plot(torque[0:38,0],torque[0:38,1])
pl.plot(torque[39:39*2,0],torque[39:39*2,1])
pl.title('torque')
pl.show()
