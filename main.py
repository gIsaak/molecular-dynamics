import numpy as np
import matplotlib.pyplot as plt
import time

##### Parameters #####

# Lennard-Jones parameters Argon
#eps = 119.8 # k_b
#sigma = 3.405 # Angstrom
#mass = 39.948*1.66e-27 #Kg

L = 10 # size of each periodic unit cell L (in units of sigma)

numOfParticles = 2 # 2 particles
numOfDimensions = 3 # 3D

num_t = 800 # time steps
timestep = 0.001 # time between iterations of Euler's method


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

################# Begin main program ########################

# set random starting point for particles
# Particle 1
setPXcoord(L/2,0,0)
setPYcoord(L/3,0,0)
setPZcoord(0,0,0)
setPXvel(10,0,0)
setPYvel(50,0,0)
setPZvel(0,0,0)
# Particle 2
setPXcoord(1,1,0)
setPYcoord(1,1,0)
setPZcoord(0,1,0)
setPXvel(20,1,0)
setPYvel(-50,1,0)
setPZvel(0,1,0)

##### Simulation #####
fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
p1 = ax.scatter(parameterMatrix[0][0][0][0],parameterMatrix[0][1][0][0], color='r')
p2 = ax.scatter(parameterMatrix[1][0][0][0],parameterMatrix[1][1][0][0], color='b')
ax.set_xlim((0,L))
ax.set_ylim((0,L))
plt.ion()
plt.show()
plt.pause(0.001)

for j in range(num_t):
    i = j%num_t # don't go over indices of parameterMatrix
    iterateCoordinates(i)
    iterateVelocities(i)
    p1.remove()
    p2.remove()
    p1 = ax.scatter(getPXcoord(0,i),getPYcoord(0,i),color='r')
    p2 = ax.scatter(getPXcoord(1,i),getPYcoord(1,i),color='b')
    plt.pause(0.00005)
    #time.sleep(0.05)

# Plot U
#time = np.arange(0, num_t*timestep, timestep)
#plt.plot(time, U)
#plt.xlabel('time')
#plt.ylabel('U')
#plt.show()
