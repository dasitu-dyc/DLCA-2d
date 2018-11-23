# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys as sys

class Aggregate:
    def __init__(self):
        self.mParticle = 0
        self.index_list = []
        self.mass_center = np.array([0.0, 0.0])
        self.R_g = 0.0   # radius of gyration

    def move(self, dx):
        global x_coord, y_coord
        for i in self.index_list:
            x_coord[i] += dx[0]
            y_coord[i] += dx[1]

    def update_mass_center(self):
        sum_coord = np.array([0.0, 0.0])
        for i in self.index_list :
            sum_coord[0] += x_coord[i]
            sum_coord[1] += y_coord[i]
        self.mass_center = sum_coord / self.mParticle

    def update_R_g(self):
        sum_r_2 = 0.0
        for i in self.index_list :
            sum_r_2 += (x_coord[i] - self.mass_center[0])**2 + (y_coord[i] - self.mass_center[1])**2 
        self.R_g = np.sqrt(sum_r_2 / self.mParticle)

    def printing(self):
        print 'm = %d' % self.mParticle
        print 'particles:', self.index_list
        print 'mass_center:', self.mass_center
        print 'R_g = %.3f' % self.R_g
        for i in self.index_list :
            print "[%.3f, %.3f]" % (x_coord[i], y_coord[i]),
        print '\n'

def get_initial_coordinates():
    x_coord = [ np.random.random() * box_width for i in xrange(n_particles)]
    y_coord = [ np.random.random() * box_width for i in xrange(n_particles)]
    
    return x_coord, y_coord 

def get_initial_velocities():
    x_vel = [k3 * np.random.normal(0.0, 1.0) for i in xrange(n_particles)]
    y_vel = [k3 * np.random.normal(0.0, 1.0) for i in xrange(n_particles)]

    return x_vel, y_vel

def get_initial_aggregates():
    aggregate_list = [Aggregate() for i in xrange(n_particles)]
    for m in xrange(n_particles):
        aggregate_list[m].mParticle = 1
        aggregate_list[m].index_list.append(m)
        aggregate_list[m].mass_center = np.array( [x_coord[m], y_coord[m]] )
        aggregate_list[m].R_g = 0.0
    
    return aggregate_list

def update_coordinates(x_coord, y_coord):
    # update coordinates
    for m in xrange(n_particles):
        if aggregate_list[m].mParticle > 1:
            for i in aggregate_list[m].index_list :
                x_coord[i] += x_vel[i] * dt
                y_coord[i] += y_vel[i] * dt
            aggregate_list[m].update_mass_center()
        elif aggregate_list[m].mParticle == 1:
            for i in aggregate_list[m].index_list :
                x_coord[i] += x_vel[i] * dt
                y_coord[i] += y_vel[i] * dt

                if x_coord[i] < 0.0 or x_coord[i] > box_width:
                    x_coord[i] =  x_coord[i] % box_width

                if y_coord[i] < 0.0 or y_coord[i] > box_width:
                    y_coord[i] =  y_coord[i] % box_width

            aggregate_list[m].update_mass_center()

    return x_coord, y_coord

def update_velocities(x_vel, y_vel):
    for m in xrange(n_particles):
        if aggregate_list[m].mParticle > 1:
            x_vel_agg = 0.0
            y_vel_agg = 0.0
            for i in xrange(aggregate_list[m].mParticle):
                x_vel_agg += k3 * np.random.normal(0.0, 1.0)
                y_vel_agg += k3 * np.random.normal(0.0, 1.0)

            x_vel_agg = x_vel_agg / aggregate_list[m].mParticle
            y_vel_agg = y_vel_agg / aggregate_list[m].mParticle
            for j in aggregate_list[m].index_list :
                x_vel[j] = x_vel_agg
                y_vel[j] = y_vel_agg
        elif aggregate_list[m].mParticle == 1:
            for j in aggregate_list[m].index_list :
                x_vel[j] = k3 * np.random.normal(0.0, 1.0)
                y_vel[j] = k3 * np.random.normal(0.0, 1.0)

    return x_vel, y_vel

def aabb_overlay(i, j):
    if x_coord[i] - x_coord[j] >= -2.0 and x_coord[i] - x_coord[j] <= 2.0 and y_coord[i] - y_coord[j] >= -2.0 and y_coord[i] - y_coord[j] <= 2.0:
        return True
    else:
        return False

def aggregate_overlay(m, k):
    for i in aggregate_list[m].index_list:
        for j in aggregate_list[k].index_list:
                if aabb_overlay(i, j) and (x_coord[i] - x_coord[j]) ** 2 + (y_coord[i] - y_coord[j]) ** 2 < 4.0:
                    # the distance between two particles is a0, and move them backforward
                    # (x_{i,k+1}-\dot{x}_{i,k}t)^2 + (x_{j,k+1}-\dot{x}_{j,k}t)^2 = a_0^2
                    # get the positive root
                    a = (x_vel[i] - x_vel[j]) ** 2 + (y_vel[i] - y_vel[j]) ** 2
                    # ideally the b is the -2.0* inner product of relative displacment and relative velocity, and the b should be positive, but in case that the relative speed is so large that the new coordinates after dt is where the tow particles cross each other, which makes the b negetive.
                    b = - 2.0 * (x_coord[i] - x_coord[j]) * (x_vel[i] - x_vel[j]) - 2.0 * (y_coord[i] - y_coord[j]) * (y_vel[i] - y_vel[j])
                    c = (x_coord[i] - x_coord[j]) ** 2 + (y_coord[i] - y_coord[j]) ** 2 - 4.0

                    if b > 0.0: 
                        B = 4.0*a/b/b * c
                        C = -1.0 - np.sqrt(1 - B)
                        t = 2.0*c/b/C    # > 0
                    elif b< 0.0:
                        B = 4.0*a/b/b * c
                        C = -1.0 - np.sqrt(1 - B)
                        t = 0.5*b/a*C    # > 0
                    else:
                        t = 0.0

                    dxi = np.array([-x_vel[i] * t, -y_vel[i] * t])
                    dxj = np.array([-x_vel[j] * t, -y_vel[j] * t])
                    aggregate_list[m].move(dxi)
                    aggregate_list[k].move(dxj)

                    # dump kth aggregate and update mth 
                    aggregate_list[m].mParticle += aggregate_list[k].mParticle
                    aggregate_list[k].mParticle = 0
                    aggregate_list[m].index_list += aggregate_list[k].index_list
                    aggregate_list[m].update_mass_center()
                    aggregate_list[m].update_R_g()
                    global n_aggregates
                    n_aggregates -= 1
                    # if m == 19 : aggregate_list[m].printing()
                    return
                '''
                    global err
                    print m, k
                    print aggregate_list[m].index_list
                    print aggregate_list[k].index_list
                    err += 1
                    if err >5 : sys.exit(1)
                    '''

def handle_overlay():
    a = 1

def aggregation():
    for m in xrange(n_particles):
        if aggregate_list[m].mParticle > 0:
            for k in xrange(m):
                if aggregate_list[k].mParticle > 0:
                     aggregate_overlay(m, k)
                      

def test_overlay():
    x_list = []
    y_list = []
    for i in xrange(n_particles):
        for j in xrange(i):
            if aabb_overlay(i, j):
                x = (x_coord[i] + x_coord[j]) / 2
                y = (y_coord[i] + y_coord[j]) / 2
                x_list.append(x)
                y_list.append(y)

    return x_list, y_list

def revise_coord(x_coord, y_coord):
    # boundary
    x_fig=x_coord
    y_fig=y_coord
    for i in xrange(n_particles):
        if x_coord[i] < 0.0:
            x_fig[i] =  x_coord[i] % box_width
            # x_vel[i] = -x_vel[i]
        elif x_coord[i] > box_width:
            x_fig[i] =  x_coord[i] % box_width
            # x_vel[i] = -x_vel[i]

        if y_coord[i] < 0.0:
            y_fig[i] =  y_coord[i] % box_width
            #y_vel[i] = -y_vel[i]
        elif y_coord[i] > box_width:
            y_fig[i] =  y_coord[i] % box_width
            #y_vel[i] = -y_vel[i]
    return x_fig, y_fig

def refresh(framenum):
    global x_coord, y_coord
    global x_vel, y_vel
    x_coord, y_coord = update_coordinates(x_coord, y_coord)
    aggregation()
    x_list, y_list = test_overlay()
    x_vel, y_vel = update_velocities(x_vel, y_vel)
    x_fig, y_fig = revise_coord(x_coord, y_coord)
    plt.cla()
    coord = plt.plot(x_fig, y_fig, 'ro')
    overlay = plt.plot(x_list, y_list, 'b*')
#    vel = plt.quiver(x_coord, y_coord, x_vel, y_vel)
    plt.plot([0, box_width, box_width, 0, 0], [0, 0, box_width, box_width, 0])

    return framenum


np.random.seed(20181122)
n_particles = 30
box_width = 50
n_steps = 100

# characteristic length
a0 = 10e-9    # radius of primary particle (10nm): m
k_B = 1.38064853e-23    # Boltzmann constant: J*K^-1
T = 300    # absolute temperature: K
mu_D_water = 8.9e-4    # dynamic viscosity of water at 25C: Pa*s
gamma = 6 * np.pi * mu_D_water * a0    # friction coefficient: N*m^-1*s
gauss_k = np.sqrt(2 * k_B / gamma)    # coefficient for the gauss process: m*s^-0.5

# characteristic time or a0^2/D with D the diffusion coefficient
t0 = gamma * a0 * a0 / k_B / T
dt = 0.1 * t0
k3 = np.sqrt(2 * k_B * T / a0 / a0 / gamma/ dt)

global x_coord, y_coord
global x_vel, y_vel
global err
err = 0
global n_aggregates
n_aggregates = n_particles

x_coord, y_coord = get_initial_coordinates()
x_vel, y_vel = get_initial_velocities()

aggregate_list = get_initial_aggregates()

'''
while n_aggregates > 100:
    x_coord, y_coord = update_coordinates(x_coord, y_coord)
    aggregation()
    x_vel, y_vel = update_velocities(x_vel, y_vel)

for m in xrange(n_particles):
    if aggregate_list[m].mParticle > 50 :
        print "%f %d" % (aggregate_list[m].R_g, aggregate_list[m].mParticle)
'''
fig, ax = plt.subplots()
## plot a final pictrue
# coord = plt.plot(x_coord, y_coord, 'ro')
# vel = plt.quiver(x_coord, y_coord, x_vel, y_vel)
# plt.plot([0, box_width, box_width, 0, 0], [0, 0, box_width, box_width, 0])

# ani = animation.FuncAnimation(fig, refresh, frames=100, interval=10, repeat=True)    # a repeat anime
ani = animation.FuncAnimation(fig, refresh, frames=n_steps, interval=50, repeat=True)
ax.set_aspect(1.0)

plt.show()
# ani.save("particles.mp4", fps=20)
