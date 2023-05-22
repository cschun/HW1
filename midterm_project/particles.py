
import numpy as np
from numpy import random
import matplotlib.pyplot as plt

G = 4*(np.pi**2)
Rmin = 3.094e+10

def initial_conditions(M):
    x = []
    y = []
    vx = []
    vy = []
    r = []

    e = 0.2*Rmin
    
    for i in range(12):
        j = 12+(3*i)
        rad = (0.2+(i*0.5))*Rmin
        angles = np.linspace(0,1,j)*(np.pi*2)
        x.extend(rad*np.cos(angles))
        y.extend(rad*np.sin(angles))
        r.extend(np.full(shape=j, fill_value=rad))
        
        v0 = np.sqrt(G*M*rad/(rad**2+e**2))
        vx.extend(v0*np.cos(angles))
        vy.extend(v0*np.sin(angles))

    r = np.array(r, float)

    #plt.scatter(x,y)
    #plt.savefig("initial.jpg")

    return x, y, vx, vy, r
    sourceFile = open('initial.txt', 'w')
    for i in range(len(r)):
        print('{:f}\t{:f}\t{:f}\t{:f}'.format(x[i], y[i], v[i], r[i]), file = sourceFile)
    sourceFile.close()

def acceleration(pos, M1, M2):
    x = [0, 10, pos[0]]
    y = [0, 10, pos[1]]
    z = [0, 0, pos[2]]
    
    ax=0
    ay=0
    az=0

    m = [M1, M2, 0]

    for i in range(2):
        dx = x[i] - x[2]
        dy = y[i] - y[2]
        dz = z[i] - z[2]
        r = np.sqrt(dx**2+dy**2+dz**2+(0.2*Rmin)**2)
        
        ax += -G*m[i]*dx/(r**3)
        ay += -G*m[i]*dy/(r**3)
        az += -G*m[i]*dz/(r**3)
    return np.array([ax, ay, az])

def threebody(pos, v, dt, M1, M2):
    a = acceleration(pos, M1, M2)
    v =  v + 0.5*dt*a
    pos = pos + dt*v
    a = acceleration(pos, M1, M2)
    v = v + 0.5*dt*a

    return pos, v

if __name__ == "__main__":
    M1 = 10e8
    M2 = 10e11
    x, y, vx, vy, r = initial_conditions(M1)
    z = 0

    dt = 100
    N = 101

    for j in range(1):
        pos_save = np.zeros((3,3, N+1))
        vel_save = np.zeros((3,3, N+1))

        pos_save[0] = [x[0], y[0], z]
        vel_save[0] = [vx[0], vy[0], z]
        
        for i in range(N):
            a = acceleration(pos, M1, M2)
            v =  v + 0.5*dt*a
            pos = pos + dt*v
            a = acceleration(pos, M1, M2)
            v = v + 0.5*dt*a
            pos_save[i+1] = p
            vel_save[i+1] = l

    m = [row[0] for row in pos]
    n = [row[1] for row in pos]
    #plt.scatter(x, y)
    plt.scatter(m, n)
    
    plt.savefig("idk.jpg")
    