import numpy as np
from numpy import random
import matplotlib.pyplot as plt

e = 0.2
G = 4*(np.pi**2)
m = 1010

def accel(x):
    r = (np.power(x[0],2) + np.power(x[1],2))**0.5
    a = [-G*m*x[0] / r**3, -G*m*x[1] / r**3]
    a=np.nan_to_num(a)
    return a

def leapstep(x, v, dt, accel):
    a = accel(x)  # call acceleration code

    v = v + 0.5 * dt * a  # advance vel by half-step
    x = x + dt * v  # advance pos by full-step

    a = accel(x)  # call acceleration code

    v = v + 0.5 * dt * a  # and complete vel. step
    return x, v

x = initial_conditions()  # set initial position
v = np.zeros((2, 140))  # set initial velocity
tnow = 0.0  # set initial time

    # next, set integration parameters
mstep = 256  # number of steps to take
dt = 1.0 / 32.0  # timestep for integration

xvalue = np.zeros((mstep, 140))
yvalue = np.zeros((mstep, 140))
t = np.zeros((mstep))

    # now, loop performing integration
for nstep in range(mstep):  # loop mstep times in all
    x, v = leapstep(x, v, dt, accel)  # take integration step
    tnow = tnow + dt  # and update value of time
    xvalue[nstep]=x[0]
    yvalue[nstep]=x[1]
    print(xvalue[nstep])

plt.scatter(xvalue[25],yvalue[25])
plt.savefig("idk.jpg")
plt.close()