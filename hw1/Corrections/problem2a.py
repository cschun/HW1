"""
leapint.py: program to integrate hamiltonian system using leapfrog.
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
    m = np.zeros((MAXPNT))

    n = 10
    #sun
    p[0] = np.array([0,0])
    v[0] = np.array([0, 0])
    m[0] = 1
    
    #mercury
    p[1] = np.array([0.387,0])
    v[1] = np.array([0,9.90784])
    m[1] = 1.651e-7
    
    #venus
    p[2] = np.array([0.723, 0])
    v[2] = np.array([0,7.382396])
    m[2] = 0.000002447
    
    #earth
    p[3] = np.array([1,0])
    v[3] = np.array([0,6.32415])
    m[3] = 3.0027e-6
    
    #mars
    p[4] = np.array([1.524,0])
    v[4] = np.array([0,5.05932])
    m[4] = 3.213e-7 
    
    #jupiter
    p[5] = np.array([3.171141621592332,-3.978931489863883]) # changed jupiter and saturns initial conditions as their original initial conditions led to overlapping orbits which was incorrect.
    v[5] = np.array([2.1230289054326365594, 1.8483100835142429741])
    m[5] = 0.00004365  
    
    #saturn
    p[6] = np.array([5.585802167594440, -8.274242200189370])
    v[6] = np.array([1.5765689269250231508,1.1360851399431532993])
    m[6] = 0.0002857
    
    #uranus
    p[7] = np.array([19.229,0])
    v[7] = np.array([0,1.432])
    m[7] = 0.00004365  
    
    #neptune
    p[8] = np.array([30.058,0])
    v[8] = np.array([0,1.146])
    m[8] = 0.00005149 

    #pluto
    p[9] = np.array([39.5,0])
    v[9] = np.array([0,0.996])
    m[9] = 6.58086572e-9
    
    tnow = 0.0  # set initial time

    mstep = 150000
    nout = 100
    dt = 1/24

    for nstep in range(mstep):  # loop mstep times in all
        if nstep % nout == 0:  # if time to output state
            printstate(p, v, n, m, tnow)  # then call output routine
        leapstep(p, v, n, m, dt)  # take integration step
        tnow = tnow + dt  # and update value of time
    if mstep % nout == 0:  # if last output wanted
        printstate(p, v, n, m, tnow)  # then output last step


def leapstep(p, v, n, m, dt):
    a = np.zeros((MAXPNT,2))

    accel(a, p, m, n)

    for i in range(n):
        v[i] = v[i] + 0.5 * dt * a[i]
        
    for i in range(n):
        p[i] = p[i] + dt * v[i]

    accel(a, p, m, n)  # call acceleration code

    for i in range(n):  # loop over all points...
        v[i] = v[i] + 0.5 * dt * a[i]

def accel(a, p, m, n):
    G = 4*(np.pi**2)
    
    for i in range(n):  # loop over all points...
        if (i!=0):
            ax = 0
            ay = 0

            x = p[i][0]-p[0][0]
            y = p[i][1]-p[0][1]

            r = np.sqrt(x**2+y**2)

            ax = -G*m[0]*x/r**3
            ay = -G*m[0]*y/r**3

            a[i][0] = ax
            a[i][1] = ay

def printstate(p, v, n, m, tnow):
    for i in range(n):  # loop over all points...
        #print(p[i])
        x[i].append(p[i][0])
        y[i].append(p[i][1])
        #print(x[i])
        
        #print("%8.4f%4d%12.6f%12.6f" % (tnow, i, x[i], v[i]))
    graphstate(x, y,vx, vy, tnow)
       
def graphstate(x, y,vx, vy, tnow):
    fig = plt.figure(figsize=(12, 10))

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