from mpl_toolkits.mplot3d import Axes3D # important for 3d scatter plot
import numpy as np
import matplotlib.pyplot as plt
#import time

##### Parameters #####

# Lennard-Jones parameters Argon
#eps = 119.8 # k_b
#sigma = 3.405 # Angstrom
#mass = 39.948*1.66e-27 #Kg

L = 4 # size of each periodic unit cell L (in units of sigma)

numOfParticles = 2 # 2 particles
numOfDimensions = 3 # 3D

num_t = 500 # time steps
timestep = 0.005 # time between iterations of Euler's and Verlet's method


# Create n x d x 3 numpy array of floats "PC3T" to store n particles
# P = particle, C = coordinate (x, y or z), 3 = (position,velocity, force), T = timestep
# in d dimensions with 3 (position, velocity and force) parameters, and num_t timesteps stored.

# For example, access particle 2's position in the y-direction at time_step=0 as
# PC3T[1][1][0][0]
# Access particle 1's momentum in the x-direction at time_step=1 as
# PC3T[0][0][1][1]

PC3T = np.zeros((numOfParticles,numOfDimensions,3,num_t), dtype=float)


# Initialize potential energy matrix U and total energy E
# TODO maybe make matrix?

U = np.zeros((num_t,1), dtype=float)
E = np.zeros((num_t,1), dtype=float)

# Get particles positions and velocities
# particle index p in 0,..,n-1
def getPXcoord(p,ts):
    return PC3T[p][0][0][ts]
def getPYcoord(p,ts):
    return PC3T[p][1][0][ts]
def getPZcoord(p,ts):
    return PC3T[p][2][0][ts]
def getPXvel(p,ts):
    return PC3T[p][0][1][ts]
def getPYvel(p,ts):
    return PC3T[p][1][1][ts]
def getPZvel(p,ts):
    return PC3T[p][2][1][ts]

# Set particles positions and velocities
# particle index p in 0,..,n-1
def setPXcoord(val,p,ts):
    PC3T[p][0][0][ts] = val
def setPYcoord(val,p,ts):
    PC3T[p][1][0][ts] = val
def setPZcoord(val,p,ts):
    PC3T[p][2][0][ts] = val
def setPXvel(val,p,ts):
    PC3T[p][0][1][ts] = val
def setPYvel(val,p,ts):
    PC3T[p][1][1][ts] = val
def setPZvel(val,p,ts):
    PC3T[p][2][1][ts] = val


def getParticleDistance(p1,p2,ts):
    # function that returns distance between p1 and p2
    # at the given timestep "ts"
    # (returned as a 4-tuple: r,x1-2,y1-2,z1-2)
    # includes periodic boundary conditions
    x_dist = getPXcoord(p1,ts) - getPXcoord(p2,ts)
    y_dist = getPYcoord(p1,ts) - getPYcoord(p2,ts)
    z_dist = getPZcoord(p1,ts) - getPZcoord(p2,ts)
    x_dist = (x_dist + L/2)%L - L/2
    y_dist = (y_dist + L/2)%L - L/2
    z_dist = (z_dist + L/2)%L - L/2

    r = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)
    return r,x_dist,y_dist,z_dist


def getTotalEnergy(ts):
    # calculates sum of potential energy and kinetic energy at the timestep
    # of particles in a unit cell
    # Natural units
    T = 0
    for i in range(numOfParticles):
        T = T + getPXvel(i,ts)**2 + getPYvel(i,ts)**2 + getPZvel(i,ts)**2
    E[ts] = U[ts] + T/2


def getForce(ts):
    # calculates the force on each particle and
    # U of the system  at a given timestep

    for i in range(numOfParticles):
        j = 0
        while j < i:
            r,rel_x,rel_y,rel_z = getParticleDistance(i,j,ts)
            invr6 = (1/r)**6 #precomputes (1/r)**6
            grad = 24/r * (-2*invr6**2  + invr6)
            # Compute forces
            PC3T[i][0][2][ts] = PC3T[i][0][2][ts] - grad*rel_x/r  #fx particle i
            PC3T[i][1][2][ts] = PC3T[i][1][2][ts] - grad*rel_y/r  #fy particle i
            PC3T[i][2][2][ts] = PC3T[i][2][2][ts] - grad*rel_z/r  #fz particle i
            PC3T[j][0][2][ts] = PC3T[j][0][2][ts] + grad*rel_x/r  #fx particle j
            PC3T[j][1][2][ts] = PC3T[j][1][2][ts] + grad*rel_y/r  #fy particle j
            PC3T[j][2][2][ts] = PC3T[j][2][2][ts] + grad*rel_z/r  #fz particle j
            j += 1
            # Compute U
            U[ts] = 4*(invr6**2 - invr6)

# Euler algorithm
def iterateCoordinates(ts):
    # takes current position and velocity at timestep ts and updates particle coordinates for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    for i in range(numOfParticles):
        newPXcoord = (getPXcoord(i,ts) + getPXvel(i,ts)*timestep)%L
        newPYcoord = (getPYcoord(i,ts) + getPYvel(i,ts)*timestep)%L
        newPZcoord = (getPZcoord(i,ts) + getPZvel(i,ts)*timestep)%L

        setPXcoord(newPXcoord,i,next_ts)
        setPYcoord(newPYcoord,i,next_ts)
        setPZcoord(newPZcoord,i,next_ts)


def iterateVelocities(ts):
    # takes current velocity and force at timestep ts and updates particle velicities for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    getForce(ts)
    for i in range(numOfParticles):
        newPXvel = getPXvel(i,ts) + PC3T[i][0][2][ts]*timestep
        newPYvel = getPYvel(i,ts) + PC3T[i][1][2][ts]*timestep
        newPZvel = getPZvel(i,ts) + PC3T[i][2][2][ts]*timestep

        setPXvel(newPXvel,i,next_ts)
        setPYvel(newPYvel,i,next_ts)
        setPZvel(newPZvel,i,next_ts)


# Velocity-Verlet algorithm
def iterateCoordinates_Verlet(ts):
    # takes current position, velocity and force at timestep ts and
    # updates particle coordinates for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1
    # To avoid redundant computations getForce is only called once per timestep in
    # iterateVelocities_Verlet
    #getForce(ts)
    for i in range(numOfParticles):
        newPXcoord = (getPXcoord(i,ts) + getPXvel(i,ts)*timestep + \
                            0.5*PC3T[i][0][2][ts]*timestep**2)%L
        newPYcoord = (getPYcoord(i,ts) + getPYvel(i,ts)*timestep + \
                            0.5*PC3T[i][1][2][ts]*timestep**2)%L
        newPZcoord = (getPZcoord(i,ts) + getPZvel(i,ts)*timestep + \
                            0.5*PC3T[i][2][2][ts]*timestep**2)%L

        setPXcoord(newPXcoord,i,next_ts)
        setPYcoord(newPYcoord,i,next_ts)
        setPZcoord(newPZcoord,i,next_ts)

def iterateVelocities_Verlet(ts):
    # takes current velocity and force at timestep ts and ts+1
    # and updates particle velicities for timestep ts+1
    # changed if condition to ts+2 due to need to force at ts+1 for updating velocities
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    # Position at time ts+1 should already be stored in memory by previous call of
    # iterateCoordinates_Verlet(ts)

    # Get force at time ts+1
    getForce(ts+1)
    for i in range(numOfParticles):
        newPXvel = getPXvel(i,ts) + 0.5*timestep*(PC3T[i][0][2][ts] + PC3T[i][0][2][ts+1])
        newPYvel = getPYvel(i,ts) + 0.5*timestep*(PC3T[i][1][2][ts] + PC3T[i][1][2][ts+1])
        newPZvel = getPZvel(i,ts) + 0.5*timestep*(PC3T[i][2][2][ts] + PC3T[i][2][2][ts+1])

        setPXvel(newPXvel,i,next_ts)
        setPYvel(newPYvel,i,next_ts)
        setPZvel(newPZvel,i,next_ts)

def initializeParticles(ts):
    # creates a random position and velocity for each particle at the given
    # timestep ts
    #
    # particle positions are generated to be distributed around the box such that no 2 particles
    # begin to close to each other; this is done by slicing the box in the first
    # dimension (L/numOfParticles) and only allowing 1 particle to exist in each slice
    #
    # particle velocities are chosen randomly in magnitude and direction, with
    # the condition that no initial velocity is faster than some fraction the width of
    # the box per timestep
    
    # this variable limits how fast the particles can initially move
    # for example, 10 means that the particles can move no further than 1/10
    # the length of the box per timestep
    maxInitialSpeedParameter = 20
    
    dim_components = np.zeros((numOfDimensions,1),dtype=float)
    
    for i in range(numOfParticles):
        # first generate positions in each dimension
        dim_components = np.random.rand(numOfDimensions,1)*L # scale to sizeOfBox
        dim_components[0] = dim_components[0]/numOfParticles/1.1 + i/L # slice in dimension 1 to separate particles initially (1.1 gives space between "sliced" regions)
        # pass vector to method to fill parameter matrix:
        # arguments: random vector, particle number, coord/vel, timestep
        addToParameterMatrix(dim_components,i,0,ts)
        
        # next generate velocities in each dimension, limited according to maxInitialSpeedParameter above
        dim_components = np.random.rand(numOfDimensions,1)/np.sqrt(numOfDimensions)*L/timestep/maxInitialSpeedParameter
        addToParameterMatrix(dim_components,i,1,ts)

def addToParameterMatrix(dim_components,pnum,xv,ts):
    # treat this  as Protected method
    #
    # function called by initializePartciels() to load randomly generated initial
    # positions/velocities into parameter matrix
    
    # TODO
    pass


################# Begin main program ########################

##### Set initial positions/velocities for all particles ####

initializeParticles(0)

# Particle 1
setPXcoord(L/2,0,0)
setPYcoord(L/3,0,0)
setPZcoord(0,0,0)
setPXvel(1,0,0)
setPYvel(2.4,0,0)
setPZvel(0,0,0)
# Particle 2
setPXcoord(1,1,0)
setPYcoord(1,1,0)
setPZcoord(0,1,0)
setPXvel(2,1,0)
setPYvel(-5,1,0)
setPZvel(0,1,0)


##### Simulation #####
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p1 = ax.scatter(getPXcoord(0,0),getPYcoord(0,0),getPZcoord(0,0), color='r')
p2 = ax.scatter(getPXcoord(1,0),getPYcoord(1,0),getPZcoord(1,0), color='b')
ax.set_xlim((0,L))
ax.set_ylim((0,L))
ax.set_zlim((0,L))
plt.ion()
plt.show()
plt.pause(0.01)

# vp1 and vp2 to keep track of the velocities
vp1 = np.zeros(num_t)
vp2 = np.zeros(num_t)

# keep track of particle distance, potential energy, and kinetic energy
particleDistances = np.zeros(num_t)

# TODO figure out a solution for the wrap around of timestep index

# use iteration in range(num_t) and iterateCoordinates for euler method
# use iteration in range(num_t-1) due to need of "future" force in Verlet algorithm

# TODO find a better implementation of this

# REMOVE THIS FOR EULER METHOD
# Get forces for initial positions
getForce(0)

for j in range(num_t-1):
    i = j%(num_t-1) # don't go over indices of PC3

    #Euler
    #iterateCoordinates(i)
    #iterateVelocities(i)

    #Velocity-Verlet
    iterateCoordinates_Verlet(i)
    iterateVelocities_Verlet(i)

    p1.remove()
    p2.remove()
    p1 = ax.scatter(getPXcoord(0,i),getPYcoord(0,i), getPZcoord(0,i),color='r')
    p2 = ax.scatter(getPXcoord(1,i),getPYcoord(1,i), getPZcoord(0,i),color='b')
    plt.pause(0.000005)

    vp1[i] = getPXvel(0,i)
    vp2[i] = getPXvel(1,i)
    particleDistances[i],a,b,c = getParticleDistance(0,1,i)
    #time.sleep(0.05)


##### Plots #########
plt.ioff()
time = np.arange(0, num_t*timestep, timestep)
plot_fig,a = plt.subplots(2,2)
a[0][0].plot(time,particleDistances,color='r',label = 'Inter-particle distance')
a[0][0].set_ylabel('Particle Distance')
a[0][0].set_xlabel('time')
a[0][1].plot(time, U,color='m',label='Potential Energy')
a[0][1].set_ylabel('Potential Energy')
a[0][1].set_xlabel('time')
a[1][0].plot(time, E-U,color='g',label='Kinetic Energy')
a[1][0].set_ylabel('Kinetic Energy')
a[1][0].set_xlabel('time')
a[1][1].plot(time, E,color='b',label='Total Energy')
a[1][1].set_ylabel('Total Energy')
a[1][1].set_xlabel('time')
plt.show()
