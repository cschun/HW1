import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


t = 0  # current time of the simulation
tEnd = .5   # time at which simulation ends
dt = 0.001  # timestep
#rmin in kpc
rmin = 25
softening = 0.2*rmin  # softening length
m1 = 10e11
m2 = 10e11
#G in terms of solar masses and kpc 
G = 4.300917270038e-06
e = 0.5

#Functions to get initial values in r,phi
def r(phi):
    return rmin * (1 + e) / (1 - e * np.cos(phi))

def vr(phi):
    return -e * np.sqrt(G * (m1 + m2) / (rmin * (1 + e))) * np.sin(phi)

def vphi(phi):
    return np.sqrt(G * (m1 + m2) / (rmin * (1 + e))) * (1 - e * np.cos(phi))
#Converting to x and y 
def vx(phi):
    return vr(phi) * np.cos(phi) - vphi(phi) * np.sin(phi)

def vy(phi):
    return vr(phi) * np.sin(phi) + vphi(phi) * np.cos(phi)

def x(phi):
    return r(phi) * np.cos(phi)

def y(phi):
    return r(phi) * np.sin(phi)
#Obtain specific x and y for each mass
def x1(phi):
    return m2 / (m1 + m2) * x(phi)

def y1(phi):
    return m2 / (m1 + m2) * y(phi)

def x2(phi):
    return -m1 / (m1 + m2) * x(phi)

def y2(phi):
    return -m1 / (m1 + m2) * y(phi)
#Obtain specific velocities for each mass
def vx1(phi):
    return m2 / (m1 + m2) * vx(phi)

def vy1(phi):
    return m2 / (m1 + m2) * vy(phi)

def vx2(phi):
    return -m1 / (m1 + m2) * vx(phi)

def vy2(phi):
    return -m1 / (m1 + m2) * vy(phi)
#Using the provided acceleration code for leapfrog
def getAcc(pos, mass, G, softening):
    """
    Calculate the acceleration on each particle due to Newton's Law
        pos  is an N x 3 matrix of positions
        mass is an N x 1 vector of masses
        G is Newton's Gravitational constant
        softening is the softening length
        a is N x 3 matrix of accelerations
    """
    # positions r = [x,y,z] for all particles
    x = pos[:, 0:1]
    y = pos[:, 1:2]
    z = pos[:, 2:3]

    # matrix that stores all pairwise particle separations: r_j - r_i
    dx = x.T - x
    dy = y.T - y
    dz = z.T - z

    # matrix that stores 1/r^3 for all particle pairwise particle separations
    inv_r3 = dx**2 + dy**2 + dz**2 + softening**2
    inv_r3[inv_r3 > 0] = inv_r3[inv_r3 > 0] ** (-1.5)

    ax = G * (dx * inv_r3) @ mass
    ay = G * (dy * inv_r3) @ mass
    az = G * (dz * inv_r3) @ mass

    # pack together the acceleration components
    a = np.hstack((ax, ay, az))

    return a

def leapfrog(acc, vel, pos, Nt, t):
    for i in range(Nt):
        # (1/2) kick
        vel += acc * dt / 2.0

        # drift
        pos += vel * dt

        # update accelerations
        acc = getAcc(pos, mass, G, softening)

        # (1/2) kick
        vel += acc * dt / 2.0

        # update time
        t += dt

        # save energies, positions for plotting trail
        pos_save[:, :, i + 1] = pos
    return pos_save

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

def animation_graph(xx, yy):
    fig = plt.figure()
    ax = plt.axes(xlim=(-40, 40), ylim=(-30, 30), xlabel="x", ylabel="v")
    line, = ax.plot([], [], lw=2, marker=".")
    line2, = ax.plot([], [], lw=2, marker=".")
    line3, = ax.plot([], [], lw=2, marker=".")
    line4, = ax.plot([], [], lw=2, marker=".")

    def init():
        line.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        line4.set_data([], [])
        return line, line2, line3, line4
    def animate(i):
        line.set_data(xx[0][i], yy[0][i])
        line2.set_data(xx[1][i], yy[1][i])
        line3.set_data(xx[2][i], yy[2][i])
        line4.set_data(xx[3][i], yy[3][i])
        return line, line2, line3, line4
    
    anim = FuncAnimation(fig, animate, init_func = init, frames=len(xx[0]), interval=100)
    anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plt.close()

def static_graph(xx, yy):
    plt.scatter(xx[1], yy[1],  s=5)
    plt.scatter(xx[2], yy[2],  s=5)
    plt.savefig("2body_orbit.png")

def rotation_matrix(i):
    a = 1-np.cos(-i)
    c = np.cos(-i)
    s = np.cos(-i)
    x = np.cos(-90)
    y = np.sin(-90)
    z = 0
    rotation = [[a*x**2+c, a*x*y-s*z, a*x*z+s*y],
                [a*x*y+s*z, a*y**2+c, a*z*y-s*x],
                [a*x*z*-s*y, a*y*z+s*x, a*z**2+c]]
    return rotation

if __name__ == "__main__":
    N = 3
    np.random.seed(17)  # set the random number generator seed

    m1pos, m1vel, m1r = initial_conditions(m1)
    m2pos, m2vel, m2r = initial_conditions(m2)

    m1pos = [x - [x1(0), y1(0), 0] for x in m1pos]
    m2pos = [x - [x2(0), y2(0), 0] for x in m1pos]

    m1vel = [x + [vx1(0), vy1(0), 0] for x in m1vel]
    m2vel = [x + [vx2(0), vy2(0), 0] for x in m2vel]

    #setting up conditions
    Nt = int(np.ceil(tEnd / dt))

    pos_array = np.zeros((len(m1pos)+2, 3, Nt+1))

    for i in range(2):
        mass = np.array([[m1], [m2], [0]])
        pos = np.array([[x1(0), y1(0), 0], [x2(0), y2(0), 0], m1pos[i]])
        vel = np.array([[vx1(0), vy1(0), 0], [vx2(0), vy2(0), 0], m1vel[i]])

        vel -= np.mean(mass * vel, 0) / np.mean(mass)

        acc = getAcc(pos, mass, G, softening)

        pos_save = np.zeros((N, 3, Nt + 1))
        pos_save[:, :, 0] = pos
        t_all = np.arange(Nt + 1) * dt

        pos_save = leapfrog(acc, vel, pos, Nt, t)
        pos_array[0] = pos_save[0]
        pos_array[1] = pos_save[1]
        pos_array[i+2] = pos_save[2]

    xx = pos_array[:, 0, :]
    yy = pos_array[:, 1, :]

    #animation_graph(xx, yy)
    static_graph(xx, yy)