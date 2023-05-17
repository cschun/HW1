import numpy as np

n = 9 # number of planets
x = np.zeros((3, n))
v = np.zeros((3, n))
nstep = 100000 # total number of steps
nout = 10 # intervals between outputs
# in units of AU, day, solar mass
G = 2.959159 * 10**(-4) 
t = 0
dt = 1.0 # time step

f = open('initial_condition.txt')
content = f.readlines()
f.close()
for i in range(n):
    line = content[i].split()
    for j in range(3):
        x[j, i] = float(line[j])
        v[j, i] = float(line[j+3])

def accel(x):
    r = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
    return - G * 1.0 * x / r**3

def printstate(x, v, t, n):
    for i in range(n):
        print('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(i, t, x[0, i], x[1, i], x[2, i], v[0, i], v[1, i], v[2, i]))

# start leapfrogging
for i in range(nstep):
    if (i % nout == 0):
        printstate(x, v, t, n)
    # evolve v by a half step
    v += accel(x) * dt / 2
    # evolve x by a full step
    x += v * dt
    # evolve v by a half step
    v += accel(x) * dt / 2
    t += dt
if (nstep % nout == 0):
    printstate(x, v, t, n)

