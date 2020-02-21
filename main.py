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

num_t = 1000 # time steps
timestep = 0.001 # time between iterations of Euler's and Verlet's method


# Create n x d x 2 numpy array of floats "parameterMatrix" to store n particles
# in d dimensions with 2 (position and velocity) parameters, and num_t timesteps stored.

# For example, access particle 2's position in the y-direction at time_step=0 as
# parameterMatrix[1][1][0][0]
# Access particle 1's momentum in the x-direction at time_step=1 as
# parameterMatrix[0][0][1][1]

parameterMatrix = np.zeros((numOfParticles,numOfDimensions,2,num_t), dtype=float)

# Initialize potential energy matrix U and total energy E
# TODO maybe make matrix?

U = np.zeros((num_t,1), dtype=float)
E = np.zeros((num_t,1), dtype=float)

# Get particles positions and velocities
# particle index p in 0,..,n-1
def getPXcoord(p,ts):
    return parameterMatrix[p][0][0][ts]
def getPYcoord(p,ts):
    return parameterMatrix[p][1][0][ts]
def getPZcoord(p,ts):
    return parameterMatrix[p][2][0][ts]
def getPXvel(p,ts):
    return parameterMatrix[p][0][1][ts]
def getPYvel(p,ts):
    return parameterMatrix[p][1][1][ts]
def getPZvel(p,ts):
    return parameterMatrix[p][2][1][ts]

# Set particles positions and velocities
# particle index p in 0,..,n-1
def setPXcoord(val,p,ts):
    parameterMatrix[p][0][0][ts] = val
def setPYcoord(val,p,ts):
    parameterMatrix[p][1][0][ts] = val
def setPZcoord(val,p,ts):
    parameterMatrix[p][2][0][ts] = val
def setPXvel(val,p,ts):
    parameterMatrix[p][0][1][ts] = val
def setPYvel(val,p,ts):
    parameterMatrix[p][1][1][ts] = val
def setPZvel(val,p,ts):
    parameterMatrix[p][2][1][ts] = val

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

# TODO choose a better name
def getForce(ts):
    # calculates the force on each particle and
    # U of the system  at a given timestep
    # returns n x d numpy array of floats "forceMatrix" to store fx,fy,fz
    # for each particle
    forceMatrix = np.zeros((numOfParticles,numOfDimensions), dtype=float)

    for i in range(numOfParticles):
        j = 0
        while j < i:
            r,rel_x,rel_y,rel_z = getParticleDistance(i,j,ts)
            invr6 = (1/r)**6 #precomputes (1/r)**6
            grad = 24/r * (-2*invr6**2  + invr6)
            # Compute forces
            forceMatrix[i][0] = forceMatrix[i][0] - grad*rel_x/r
            forceMatrix[i][1] = forceMatrix[i][1] - grad*rel_y/r
            forceMatrix[i][2] = forceMatrix[i][2] - grad*rel_z/r
            forceMatrix[j][0] = forceMatrix[j][0] + grad*rel_x/r
            forceMatrix[j][1] = forceMatrix[j][1] + grad*rel_y/r
            forceMatrix[j][2] = forceMatrix[j][2] + grad*rel_z/r
            # Compute U
            U[ts] = 4*(invr6**2 - invr6)

            j += 1
    return forceMatrix


def iterateCoordinates(ts):
    # takes current position and velocity at timestep ts and updates particle coordinates for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    for i in range(numOfParticles):
        newPXcoord = fixPos(getPXcoord(i,ts) + getPXvel(i,ts)*timestep)
        newPYcoord = fixPos(getPYcoord(i,ts) + getPYvel(i,ts)*timestep)
        newPZcoord = fixPos(getPZcoord(i,ts) + getPZvel(i,ts)*timestep)

        setPXcoord(newPXcoord,i,next_ts)
        setPYcoord(newPYcoord,i,next_ts)
        setPZcoord(newPZcoord,i,next_ts)

def fixPos(val):
    # accepts a coordinate and makes sure the value is within the interval [0,1)
    return val%L

def iterateVelocities(ts):
    # takes current velocity and force at timestep ts and updates particle velicities for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    force = getForce(ts)
    for i in range(numOfParticles):
        newPXvel = getPXvel(i,ts) + force[i][0]*timestep
        newPYvel = getPYvel(i,ts) + force[i][1]*timestep
        newPZvel = getPZvel(i,ts) + force[i][2]*timestep

        setPXvel(newPXvel,i,next_ts)
        setPYvel(newPYvel,i,next_ts)
        setPZvel(newPZvel,i,next_ts)
        
# using the same iterator structure but velocity-verlet algorithm
# ==================================
# TO DO:      
# ================================== 
# take future force into account. Now when clculating the force at one time step
# ahead we end up using the present particles position instead of the future one
# need to implement a way that calculates the position one step ahead seperately
# could use conditional function, open to discuss
        
def iterateCoordinates_Verlet(ts):
    # takes current position, velocity and force at timestep ts and 
    # updates particle coordinates for timestep ts+1
    if ts + 1 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1
        
    force = getForce(ts)
    for i in range(numOfParticles):
        newPXcoord = fixPos(getPXcoord(i,ts) + getPXvel(i,ts)*timestep + \
                            0.5*force[i][0]*timestep**2)
        newPYcoord = fixPos(getPYcoord(i,ts) + getPYvel(i,ts)*timestep + \
                            0.5*force[i][1]*timestep**2)
        newPZcoord = fixPos(getPZcoord(i,ts) + getPZvel(i,ts)*timestep + \
                            0.5*force[i][0]*timestep**2)

        setPXcoord(newPXcoord,i,next_ts)
        setPYcoord(newPYcoord,i,next_ts)
        setPZcoord(newPZcoord,i,next_ts)
        
def iterateVelocities_Verlet(ts):
    # takes current velocity and force at timestep ts and ts+1
    # and updates particle velicities for timestep ts+1
    # changed if condition to ts+2 due to need to force at ts+1 for updating velocities
    if ts + 2 == num_t:
        next_ts = 0
    else:
        next_ts = ts + 1

    force = getForce(ts)
    # iterateCoordinates_Verlet(ts)
    force_1 = getForce(ts+1)
    for i in range(numOfParticles):
        newPXvel = getPXvel(i,ts) + 0.5*timestep*(force_1[i][0] + force[i][0])
        newPYvel = getPYvel(i,ts) + 0.5*timestep*(force_1[i][1] + force[i][0])
        newPZvel = getPZvel(i,ts) + 0.5*timestep*(force_1[i][2] + force[i][0])

        setPXvel(newPXvel,i,next_ts)
        setPYvel(newPYvel,i,next_ts)
        setPZvel(newPZvel,i,next_ts)

################# Begin main program ########################

# set random starting point for particles
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

# use iteration in range(num_t) and iterateCoordinates for euler method

# use iteration in range(num_t-1) due to need of "future" force in Verlet 
# and iterateCoordinates_Verlet; 
for j in range(num_t):
    i = j%num_t # don't go over indices of parameterMatrix
    iterateCoordinates(i)
    iterateVelocities(i)
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
# Distance plot
time = np.arange(0, num_t*timestep, timestep)
plot_fig,a = plt.subplots(4,1)
a[0].plot(time,particleDistances,color='r',label = 'Inter-particle distance')
a[0].set_ylabel('Particle Distance')
a[1].plot(time, U,color='m',label='Potential Energy')
a[1].set_ylabel('Potential Energy')
a[2].plot(time, E-U,color='g',label='Kinetic Energy')
a[2].set_ylabel('Kinetic Energy')
a[3].plot(time, E,color='b',label='Total Energy')
a[3].set_ylabel('Total Energy')

plt.xlabel('time')