"""
problem2a.py. simulation of planets orbitting the sun
the unit system used in this is AU-Years-Solar Mass
"""
import numpy as np
import matplotlib.pyplot as plt

MAXPNT = 10
x = [[],[],[],[],[],[],[],[],[],[]]
y = [[],[],[],[],[],[],[],[],[],[]]
vx = [[],[],[],[],[],[],[],[],[],[]]
vy = [[],[],[],[],[],[],[],[],[],[]]
    
def main():
    p = np.zeros((MAXPNT,2))
    v = np.zeros((MAXPNT,2))
    n = 10
    #sun
    p[0] = np.array([0,0])
    v[0] = np.array([0, 0])
    
    #mercury
    p[1] = np.array([0.387,0])
    v[1] = np.array([0,9.90784])
    
    #venus
    p[2] = np.array([0.723, 0])
    v[2] = np.array([0,7.382396])
    
    #earth
    p[3] = np.array([1,0])
    v[3] = np.array([0,6.32415])
    
    #mars
    p[4] = np.array([1.524,0])
    v[4] = np.array([0,5.05932])
    
    #jupiter
    p[5] = np.array([5.2,0])
    v[5] = np.array([0,2.7552232])
    
    #saturn
    p[6] = np.array([9.538,0])
    v[6] = np.array([0,1.432])
    
    #uranus
    p[7] = np.array([19.229,0])
    v[7] = np.array([0,1.432])
    
    #neptune
    p[8] = np.array([30.058,0])
    v[8] = np.array([0,1.146])

    #pluto
    p[9] = np.array([39.5,0])
    v[9] = np.array([0,0.996])
    
    tnow = 0.0  # set initial time

    mstep = 150000
    nout = 100
    dt = 1/24

    for nstep in range(mstep):  # loop mstep times in all
        if nstep % nout == 0:  # if time to output state
            printstate(p, n, tnow)  # then call output routine
        leapstep(p, v, n, dt)  # take integration step
        tnow = tnow + dt  # and update value of time
    if mstep % nout == 0:  # if last output wanted
        printstate(p, n, tnow)  # then output last step


def leapstep(p, v, n, dt):
    a = np.zeros((MAXPNT,2))

    accel(a, p, n)

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

def accel(a, p, n):
    G = 4*(np.pi**2)
    
    for i in range(n):  # loop over all points...
        if (i!=0):
            ax = 0
            ay = 0

            x = p[i][0]-p[0][0]
            y = p[i][1]-p[0][1]
            v=[x,y]
            
            r = np.sqrt(x**2+y**2)

            ax = -G*x/r**3
            ay = -G*y/r**3

            a[i][0] = ax
            a[i][1] = ay

def printstate(p, n, tnow):
    for i in range(n):  # loop over all points...
        x[i].append(p[i][0])
        y[i].append(p[i][1])
        
    graphstate(x, y, tnow)
       
def graphstate(x, y, tnow):
    time = np.linspace(0, tnow, len(x))

    plt.plot(x[0::10], y[0::10], "o", linestyle = "--", color ='yellow', markersize=10)
    plt.plot(x[1::10], y[1::10], "o", linestyle = "--", color ='#1a1a1a', markersize=1)
    plt.plot(x[2::10], y[2::10], "o", linestyle = "--", color ='#e6e6e6', markersize=1)
    plt.plot(x[3::10], y[3::10], "o", linestyle = "--", color ='#2f6a69', markersize=1)
    plt.plot(x[4::10], y[4::10], "o", linestyle = "--", color ='#993d00', markersize=1)
    plt.plot(x[5::10], y[5::10], "o", linestyle = "--", color ='#b07f35', markersize=4)
    plt.plot(x[6::10], y[6::10], "o", linestyle = "--", color ='#b08f36', markersize=4)
    plt.plot(x[7::10], y[7::10], "o", linestyle = "--", color ='#5580aa', markersize=4)
    plt.plot(x[8::10], y[8::10], "o", linestyle = "--", color ='#366896', markersize=4)
    plt.plot(x[9::10], y[9::10], "o", linestyle = "--", color ='#9ca6b7', markersize=4)
    
    plt.show()
    plt.savefig('2a.png')
    plt.close()

if __name__ == "__main__":
    main()