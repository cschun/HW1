import numpy as np

N = 2
tEnd = 1   # time at which simulation ends
dt = 1  # timestep
rmin = 25
G = 4.300917270038e-06
softening = 0.2*rmin 

def initial_conditions(M):
    pos = []
    vel = []
    r = []
    
    for i in range(12):
        j = 12+(3*i)
        rad = (0.2+(i*0.5))*rmin
        angles = np.linspace(0,1,j)*(np.pi*2)

        pos.extend(np.stack((rad*np.cos(angles), rad*np.sin(angles), np.zeros(len(angles))), axis = 1))
        r.extend(np.full(shape=j, fill_value=rad))
        
        v0 = np.sqrt(G*M*rad/(rad**2+softening**2))
        vel.extend(np.stack((v0*np.cos(angles), v0*np.sin(angles), np.zeros(len(angles))), axis = 1))


    r = np.array(r, float)

    return pos, vel, r

pos, vel, r = initial_conditions(1e11)
print(pos[0])
