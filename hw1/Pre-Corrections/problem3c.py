"""
leapint.py: program to integrate hamiltonian system using leapfrog.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython import display

MAXPNT = 10

ex = [[],[],[]]
ey = [[],[],[]]
evx = [[],[],[]]
evy = [[],[],[]]
rx = [[],[],[]]
ry = [[],[],[]]

def main():
    p = np.zeros((MAXPNT,2))
    v = np.zeros((MAXPNT,2))
    e = np.zeros((MAXPNT, 2))
    m = np.zeros((MAXPNT))

    n = 3
    #sun
    p[0] = np.array([0,0])
    v[0] = np.array([0, 0])
    e[0] = np.array([0, 0])
    m[0] = 1
    
    #earth
    p[1] = np.array([1,0])
    v[1] = np.array([0,6.32415])
    e[1] = np.array([0.01671, 1])
    m[1] = 3.0027e-6
    
    #mars
    p[2] = np.array([1.524,0])
    v[2] = np.array([0,5.05932])
    e[2] = np.array([0.093, 1.524])
    m[2] = 3.213e-7 
    
    tnow = 0.0  # set initial time

    mstep = 3650
    nout = 115
    dt = 1/25

    for nstep in range(mstep):  # loop mstep times in all
        if nstep % nout == 0:  # if time to output state
            printstate(p, v, n, m, tnow)  # then call output routine
        leapstep(p, v, n, m, dt)  # take integration step
        tnow = tnow + dt  # and update value of time
    if mstep % nout == 0:  # if last output wanted
        analytic_orbit(e, n)
        printstate(p, v, n, m, tnow)  # then output last step
        graphstate(tnow)


def leapstep(p, v, n, m, dt):
    a = np.zeros((MAXPNT,2))

    accel(a, p, m, n)

    for i in range(n):
        v[i][0] = v[i][0] + 0.5 * dt * a[i][0]
        v[i][1] = v[i][1] + 0.5 * dt * a[i][1]
        
    for i in range(n):
        p[i][0] = p[i][0] + dt * v[i][0]
        p[i][1] = p[i][1] + dt * v[i][1]

    accel(a, p, m, n)  # call acceleration code

    for i in range(n):  # loop over all points...
        v[i][0] = v[i][0] + 0.5 * dt * a[i][0]
        v[i][1] = v[i][1] + 0.5 * dt * a[i][1]

def accel(a, p, m, n):
    G = 4*(np.pi**2)
    
    for i in range(n):  # loop over all points...
        if (i!=0):
            ax = 0
            ay = 0

            x = p[i][0]-p[0][0]
            y = p[i][1]-p[0][1]
            v=[x,y]
            
            v_hat = v / np.linalg.norm(v)
            r = np.sqrt(x**2+y**2)

            ax = -G*m[0]*x/r**3
            ay = -G*m[0]*y/r**3

            a[i][0] = ax
            a[i][1] = ay

def analytic_orbit(e, n):
    for i in range(n):
        for j in range(0,360, 1):
            t = j*np.pi/180
            radius = e[i][1]*(1-e[i][0]**2)/(1+(e[i][0]*np.cos(j*np.pi/180)))
            rx[i].append(radius*np.cos(t))
            ry[i].append(radius*np.sin(t))
            
def printstate(p, v, n, m, tnow):
    for i in range(n):  # loop over all points...
        ex[i].append(p[i][0])
        ey[i].append(p[i][1])
        evx[i].append(v[i][0])
        evy[i].append(v[i][1])
        
        #print("%8.4f%4d%12.6f%12.6f" % (tnow, i, ex[i], ey[i]))
       
    
def graphstate(tnow):

    fig = plt.figure()
    ax = plt.axes(xlim=(-3, 23), ylim=(-3, 3.5), xlabel="x", ylabel="y")
    line, = ax.plot([], [], lw=2, marker=".")
    line1, = ax.plot([], [], lw=2, marker=".")

    def init():
        line.set_data([], [])
        line1.set_data([], [])
        return line, line1

    def animate(i):
        line.set_data(rx[1][0:i+1], ry[1][0:i+1])
        line1.set_data(rx[0][0:i+1], ry[0][0:i+1])
        return line, line1
    
    anim = FuncAnimation(fig, animate, init_func = init, frames=len(rx[1]), interval=10)
    anim.save('3c.gif')

    plt.close()

if __name__ == "__main__":
    main()