"""
leapint.py: program to integrate hamiltonian system using leapfrog.
"""
import numpy as np
import matplotlib.pyplot as plt

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

    mstep = 365*2
    nout = 1
    dt = 1/31

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
    fig, ax = plt.subplots()
    ax.scatter(ex[0::3], ey[0::3], color ='yellow', s=5, label="sun")
    ax.scatter(ex[1::3], ey[1::3], color ='#2f6a69', s=5, label = "integrated Earth")
    #ax.scatter(ex[2::3], ey[2::3], color ='#993d00', s=5, label = "integrated Mars")
    
    ax.scatter(rx[1::3], ry[1::3], color ='green', s=5, label = "analytical Earth")
    #ax.scatter(rx[2::3], ry[2::3], color ='red', s=5, label = "analytical Mars")
    ax.legend()

    plt.savefig('2b_earth.png')
    plt.close()

if __name__ == "__main__":
    main()