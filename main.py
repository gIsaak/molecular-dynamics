import numpy as np

##### Parameters #####
numOfParticles = 2 # particle number
numOfDimensions = 2 # number of spatial dimensions
sizeOfBox = 1 # size of each periodic unit cell L (maximum allowed distance between particles)
num_t = 100 # number of timesteps to be tracked simultaneously

# Lennard-Jones parameters (all in SI)
boltzmann = 1.38e-23
epsilon = 119.8*boltzmann
sigma = 3.405e-10

# Argon parameters (SI)
mass = 39.948*1.66e-27


# Create n x d x 2 numpy array of floats "parameterMatrix" to store n particles
# in d dimensions with 2 (position and velocity) parameters, and num_t timesteps stored.

# For example, access particle 2's position in the y-direction at time_step=0 as
# parameterMatrix[1][1][0][0]
# Access particle 1's momentum in the x-direction at time_step=1 as
# parameterMatrix[0][0][1][1]

parameterMatrix = np.zeros((numOfParticles,numOfDimensions,2,num_t), dtype=float)

def getP1Xdist(ts):
    return parameterMatrix[0][0][0][ts]
def getP2Xdist(ts):
    return parameterMatrix[1][0][0][ts]
def getP1Ydist(ts):
    return parameterMatrix[0][1][0][ts]
def getP2Ydist(ts):
    return parameterMatrix[1][1][0][ts]
def getP1Xvel(ts):
    return parameterMatrix[0][0][1][ts]
def getP2Xvel(ts):
    return parameterMatrix[1][0][1][ts]
def getP1Yvel(ts):
    return parameterMatrix[0][1][1][ts]
def getP2Yvel(ts):
    return parameterMatrix[1][1][1][ts]

# since potential is only a function of distance, we can efficiently
# compute the potential for all possible distances (given periodic BCs),
# compute the 2D gradient over the space
n = 101 # number of mesh points in box
x = np.linspace(0.01,sizeOfBox,n)
xx,yy = np.meshgrid(x,x)
# xx and yy hold a 2D array of the x and y coordinates on a 2D Grid

# evaluate the Lennard Jones Potential on the grid
lj = 4*epsilon*((sigma/sqrt(xx**2+yy**2)**12-(sigma/sqrt(xx**2+yy**2))**6)


def getParticleDistance(ts):
    # function that returns the minimum distance between the two particles at the given timestep "ts"
    # includes periodic boundary conditions
    x_dist = abs(getP1Xdist(ts)-getP2Xdist(ts))
    y_dist = abs(getP1Ydist(ts)-getP2Ydist(ts))
    if x_dist > sizeOfBox/2:
        # if the x-distance is greater than half the x-length of the box, the shortest x distance wraps over the periodic BCs
        x_dist = sizeOfBox - x_dist
    if y_dist > sizeOfBox/2;
        # same reasoning as for x coordinate
        y_dist = sizeOfBox - y_dist

    # calculate and return distance shortest distance between particles
    return sqrt(x_dist**2 + y_dist**2)

def getPotentialEnergy(ts):
    # returns the potential energy (Lennard Jones) at the given timestep
    r = getParticleDistance(ts)
    return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

def getTotalEnergy(ts):
    # calculates sum of potential energy and kinetic energy at the timestep (of particles in a unit cell)
    U = getPotentialEnergy(ts)
    T = 1/2*mass*(getP1Xvel(ts)**2+getP1Yvel(ts)**2+getP2Xvel(ts)**2+getP2Yvel(ts)**2)
    return U + T

def getForce(ts):

