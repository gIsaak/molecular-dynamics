from mpl_toolkits.mplot3d import Axes3D # important for 3d scatter plot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
from itertools import combinations
##### Parameters #####

# Lennard-Jones parameters Argon
#eps = 119.8 # k_b
#sigma = 3.405 # Angstrom
#mass = 39.948*1.66e-27 #Kg


# Deprecated functions to read/write to PC3T array
'''
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
def setPncoord(val,d,p,ts):
    # set coordinate for arbitrary dimension
    PC3T[p][d][0][ts] = val
def setPXvel(val,p,ts):
    PC3T[p][0][1][ts] = val
def setPYvel(val,p,ts):
    PC3T[p][1][1][ts] = val
def setPZvel(val,p,ts):
    PC3T[p][2][1][ts] = val
def setPnvel(val,d,p,ts):
    # set velocity for arbitrary dimension
    PC3T[p][d][1][ts] = val
'''

# Deprecated function
'''
def getParticleDistance(p1,p2,ts):

    #DEPRECATED: this is done entirely in getParticleDistanceArray
    #Returns distance between p1 and p2
    #at the given timestep "ts"
    #(returned as a 4-tuple: r,x1-2,y1-2,z1-2)
    #includes periodic boundary conditions

    x_dist = getPXcoord(p1,ts) - getPXcoord(p2,ts)
    y_dist = getPYcoord(p1,ts) - getPYcoord(p2,ts)
    z_dist = getPZcoord(p1,ts) - getPZcoord(p2,ts)
    x_dist = (x_dist + L/2)%L - L/2
    y_dist = (y_dist + L/2)%L - L/2
    z_dist = (z_dist + L/2)%L - L/2

    r = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)
    return r,x_dist,y_dist,z_dist
'''


# Deprecated function
'''
def getTotalEnergy(ts):

    #DEPRECATED: always calculate total energy as sum of potential and kinetic
    #when necessary
    #Calculates sum of potential energy and kinetic energy at the timestep
    #of particles in a unit cell
    #Natural units

    T = 0
    for i in range(numOfParticles):
        T = T + getPXvel(i,ts)**2 + getPYvel(i,ts)**2 + getPZvel(i,ts)**2
    E[ts] = U[ts] + T/2
'''

def getKineticEnergy(v,numOfParticles,numOfDimensions):
    '''
    Accepts a velocity submatrix v containing all velocity components for all
    particles. The numbers of particles and dimensions are also passed.

    The funciton returns the kinetic energy as calculated from the submatrix v.
    '''
    KE = 0
    for i in range(numOfParticles):
        for j in range(numOfDimensions):
            KE = KE + 0.5*v[i,j]**2
    return KE

def getParticleDistanceArray(x,numOfParticles,numOfDimensions,boxSize):
    '''
    Accepts a position submatrix x containing all position components for all
    particles. The numbers of particles and dimensions as well as box size are
    also passed.

    The function returns an array containing the distance between each
    unique particle pair. The order in which these distances is stored is
    evident from the loops in this function. By starting at pairIndex = 0:
        for i in range(numOfParticles):
            j = 0
            while j < i:
                ...
                <use particleDistances[pairIndex,...]>
                ...
                j += 1
                pairIndex += 1

    For each row of the array particleDistances, the first column stores
    the total distance, and the subsequent columns store the distance for
    each dimension, up to a total numOfDimensions components.

    The particleDistances array is returned for use with tracking observables
    such as pair correlation and for computing the forces acting on each particle.
    '''
    # initialize array to store particle-pair distances for all unique pairings
    # the four columns contain: r (total distance), r_x, r_y, r_z, etc.
    particleDistances = np.zeros(((numOfParticles**2-numOfParticles)/2,1+numOfDimensions),dtype=float)
    pairIndex = 0

    # we will always iterate over particle pairs as follows for consistency
    for i in range(numOfParticles):
        j = 0
        while j < i:
            r_sq = 0 # total distance squared
            for k in range(numOfDimensions):
                dist_component = x[i,k] - x[j,k]
                dist_component = (dist_component + boxSize/2)%boxSize - boxSize/2
                particleDistances[pairIndex,k+1] = dist_component
                r_sq += dist_component**2
            particleDistances[pairIndex,0] = np.sqrt(r_sq)
            pairIndex += 1
            j += 1
    return particleDistances

def getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize):
    '''
    Accepts the particleDistances array computed by getParticleDistanceArray, as well
    as the numbers of particles and dimensions and box size.

    The function returns a submatrix newF containing the forces acting on the
    particles given the particleDistances array. The function also calculates
    the total potential energy in the system and returns this alongside the array.
    '''
    pairIndex = 0
    PE = 0
    newF = np.zeros((numOfParticles,numOfDimensions),dtype=float)
    for i in range(numOfParticles):
        j = 0
        while j < i:
            r = particleDistances[pairIndex,0]
            invr6 = (1/r)**6 #precomputes (1/r)**6
            grad = 24/r * (-2*invr6**2  + invr6)
            # Compute forces
            for k in range(numOfDimensions):
                newF[i,k] = newF[i,k] - grad*particleDistances[pairIndex,k+1]/r
                newF[j,k] = newF[j,k] + grad*particleDistances[pairIndex,k+1]/r

                #PC3T[i][0][2][ts] = PC3T[i][0][2][ts] - grad*rel_x/r  #fx particle i
                #PC3T[i][1][2][ts] = PC3T[i][1][2][ts] - grad*rel_y/r  #fy particle i
                #PC3T[i][2][2][ts] = PC3T[i][2][2][ts] - grad*rel_z/r  #fz particle i
                #PC3T[j][0][2][ts] = PC3T[j][0][2][ts] + grad*rel_x/r  #fx particle j
                #PC3T[j][1][2][ts] = PC3T[j][1][2][ts] + grad*rel_y/r  #fy particle j
                #PC3T[j][2][2][ts] = PC3T[j][2][2][ts] + grad*rel_z/r  #fz particle j
            pairIndex += 1
            j += 1

            # Compute Potential Energy
            PE = PE + 4*(invr6**2 - invr6)
    return newF,PE

def iterateCoordinates_Euler(x,v,numOfParticles,numOfDimensions,timestep,boxSize):
    '''
    Accepts 2 submatrices for particle positions (x) and particle velocities (v)
    as well as numbers of particles and dimensions and the size of the timestep.
    The function returns a new submatrix newX containing the Euler-time-evolved
    positions at the next timestep.
    '''

    # create array of new particle positions according to Euler time evolution
    newX = np.zeros((numOfParticles,numOfDimensions),dtype=float)
    for i in range(numOfParticles):
        for j in range(numOfDimensions):
            newX = (x[i,j] + v[i,j]*timestep)%boxSize
    return newX

def iterateVelocities_Euler(v,f,numOfParticles,numOfDimensions,timestep):
    '''
    Accepts 2 submatrices for particle velocities (v) and particle forces (f)
    as well as numbers of particles and dimensions and the size of the timestep.
    The function returns a new submatrix newV containing the Euler-time-evolved
    velocities at the next timestep.
    '''

    newV = np.zeros((numOfParticles,numOfDimensions),dtype=float)
    for i in range(numOfParticles):
        for j in range(numOfDimensions):
            newV = (v[i,j]+f[i,j]*timestep)

        #newPXvel = getPXvel(i,ts) + PC3T[i][0][2][ts]*timestep
        #newPYvel = getPYvel(i,ts) + PC3T[i][1][2][ts]*timestep
        #newPZvel = getPZvel(i,ts) + PC3T[i][2][2][ts]*timestep
    return newV

# Verlet algorithm
def iterateCoordinates_Verlet(x,v,f,numOfParticles,numOfDimensions,timestep,boxSize):
    '''
    Accepts 3 submatrices for particle positions (x), particle velocities (v), and
    particle forces (f) as well as numbers of particles and dimensions and the size of the timestep
    and the box size. The function returns a new submatrix newX containing the Verlet-time-evolved
    positions at the next timestep.
    '''
    # To avoid redundant computations getForce is only called once per timestep in
    # iterateVelocities_Verlet

    newX = np.zeros((numOfParticles,numOfDimensions),dtype=float)
    for i in range(numOfParticles):
        for j in range(numOfDimensions):
            newX[i,j] = (x[i,j] + v[i,j]*timestep + \
                            0.5*f[i,j]*timestep**2)%boxSize

            #newPXcoord = (getPXcoord(i,ts) + getPXvel(i,ts)*timestep + \
            #                    0.5*PC3T[i][0][2][ts]*timestep**2)%L
            #newPYcoord = (getPYcoord(i,ts) + getPYvel(i,ts)*timestep + \
            #                    0.5*PC3T[i][1][2][ts]*timestep**2)%L
            #newPZcoord = (getPZcoord(i,ts) + getPZvel(i,ts)*timestep + \
            #                    0.5*PC3T[i][2][2][ts]*timestep**2)%L
        #setPXcoord(newPXcoord,i,next_ts)
        #setPYcoord(newPYcoord,i,next_ts)
        #setPZcoord(newPZcoord,i,next_ts)
    return newX

# Verlet algorithm
def iterateVelocities_Verlet(v,f,nextf,numOfParticles,numOfDimensions,timestep):
    '''
    Accepts 3 submatrices for particle velocities (v), particle forces at the same
    timestep as v (f), and particle forces at the timestep after v (nextf) as well
    as numbers of particles and dimensions and the size of the timestep
    The function returns a new submatrix newV containing the Verlet-time-evolved
    velocities at the next timestep.
    '''
    newV = np.zeros((numOfParticles,numOfDimensions),dtype=float)
    for i in range(numOfParticles):
        for j in range(numOfDimensions):
            newV[i,j] = v[i,j] + 0.5*timestep*(f[i,j] + nextf[i,j])

            #newPXvel = getPXvel(i,ts) + 0.5*timestep*(PC3T[i][0][2][ts] + PC3T[i][0][2][ts+1])
            #newPYvel = getPYvel(i,ts) + 0.5*timestep*(PC3T[i][1][2][ts] + PC3T[i][1][2][ts+1])
            #newPZvel = getPZvel(i,ts) + 0.5*timestep*(PC3T[i][2][2][ts] + PC3T[i][2][2][ts+1])

            #setPXvel(newPXvel,i,next_ts)
            #setPYvel(newPYvel,i,next_ts)
            #setPZvel(newPZvel,i,next_ts)
    return newV

def init_position(a,numOfParticles,numOfDimensions):
    '''
    Accepts lattice constant (a), numOfParticles, numOfDimensions
    Returns a submatrix newX of particles positions on fcc lattice
    as well as new number of particles (newNumOfParticles) and box size (L)
    !!! ONLY WORKS IN 3D !!!

    For a given n = L/a, an fcc compatible number of particles is
    N = (n + 1)**3 + 3*(n+1)*n**2
    The initialization rounds the user specified N to the closest fcc compatible N.
    a is kept constant and L is calculated as
    L = a*n + a with updated N
    '''
    # Find a compatible n
    func = lambda x : numOfParticles - (x + 1)**3 - 3*x**3
    n = int(np.round(fsolve(func, 1.1))) #1.1 initial guess
    # Compute L and exact N
    L = a*n + a # +a to avoid putting particles on boundary
    newNumOfParticles = (n + 1)**3 + 3*(n+1)*n**2
    # Initialize PC3T (in case N is bigger than the input N)
    newX = np.zeros((newNumOfParticles,numOfDimensions), dtype=float)
    # Put particles on cubic lattice
    p = 0
    for i in range(n+1): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n+1): # loop over x
                #p = i*(n+1)**2 + j*(n+1) + k #particle index, (n+1)**3 total on the cubic lattice
                newX[p,0] = k*a + a/2  #a/2 avoids putting particles on boundary
                newX[p,1] = j*a + a/2
                newX[p,2] = i*a + a/2
                p += 1
                #setPXcoord(k*a + a/2,p,ts)
                #setPYcoord(j*a + a/2,p,ts)
                #setPZcoord(i*a + a/2,p,ts)
    # Put particles on xz faces
    for i in range(n): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n): # loop over x
                #p = (n+1)**3 + i*n*(n+1) + j*n + k #particle index, (n+1)*n**2 on x faces
                newX[p,0] = k*a + a
                newX[p,1] = j*a + a/2
                newX[p,2] = i*a + a
                p += 1
                #setPXcoord(k*a + a,p,ts)
                #setPYcoord(j*a + a/2,p,ts)
                #setPZcoord(i*a + a,p,ts)
    # Put particles on yz faces
    for i in range(n): # loop over z
        for j in range(n): # loop over y
            for k in range(n+1): # loop over x
                #p = (n+1)**3 + (n+1)*n**2 + i*n*(n+1) + j*(n+1) + k #particle index, (n+1)*n**2 on yz faces
                newX[p,0] = k*a + a/2
                newX[p,1] = j*a + a
                newX[p,2] = i*a + a
                p += 1
                #setPXcoord(k*a + a/2,p,ts)
                #setPYcoord(j*a + a,p,ts)
                #setPZcoord(i*a + a,p,ts)
    # Put particles on xy faces
    for i in range(n+1): # loop over z
        for j in range(n): # loop over y
            for k in range(n): # loop over x
                #p = (n+1)**3 + 2*(n+1)*n**2 + i*n*n + j*n + k #particle index, (n+1)*n**2 on xy faces
                newX[p,0] = k*a + a
                newX[p,1] = j*a + a
                newX[p,2] = i*a + a/2
                p += 1
                #setPXcoord(k*a + a,p,ts)
                #setPYcoord(j*a + a,p,ts)
                #setPZcoord(i*a + a/2,p,ts)
    return newX, newNumOfParticles, L

# TODO modify to work with cleaned code
def gaussVel(temp, numOfParticles, numOfDimensions):
    '''
    Accepts temperature (temp), numOfParticles, numOfDimensions
    Returns submatrix newV of velocity components per particle initialized
    according to Gaussian distribution.
    The mean value of each velocity component array is set to zero.
    '''
    newV = np.zeros((numOfParticles,numOfDimensions), dtype=float)
    # Distribution parameters
    mu, sigma = 0, np.sqrt(temp/119.8) # mean is 0 and standard deviation in Kelvin
    for i in range(numOfDimensions):
        v_i = np.random.normal(mu, sigma, numOfParticles)
        newV[:,i] = v_i - np.mean(v_i)
    #vx = np.random.normal(mu, sigma, numOfParticles)
    #vy = np.random.normal(mu, sigma, numOfParticles)
    #vz = np.random.normal(mu, sigma, numOfParticles)
    # set mean to zeros
    #vx = vx - np.mean(vx)
    #vy = vy - np.mean(vy)
    #vz = vz - np.mean(vz)
    # load into PC3T
    return newV


def getVelocityScalingFactor(KE,bathTemperature):
    '''
    Function accepts a kinetic energy and desired system temperature, and returns
    the appropriate value by which all velocity components should be scaled such
    that the system can equilibrate to the given temperature.

    The rescalingConstant should be multiplied to the present velocity components
    using np.multiply(array,rescalingConstant)
    '''
    # calculate the rescaling
    sum_of_m_vi_squared = 2*float(KE)
    rescalingConstant = np.sqrt((numOfParticles-1)*3*bathTemperature/119.8/sum_of_m_vi_squared)

    return rescalingConstant

# TODO: make in one plot
def plotEnergy(algorithm,numOfParticles,numOfTimesteps,timestep,saveFigures,U,T):
    '''
    Creates plot for Kinetic Energy, Potential energy and Total energy.
    Plot title/name/saving/etc is handled given the parameters set in main.py

    The plots will show all energies at all steps, including before the
    desired equilibrium temperature has been reached.
    '''
    # Figure name
    name = 'Eplot_N_{}p_{}'.format(numOfParticles,algorithm)
    # Three plots in horizontal
    fig, axs = plt.subplots(3, num=name)
    plt.ioff()
    time = np.arange(0, numOfTimesteps*timestep, timestep)
    # Potentil energy
    axs[0].plot(time, U,color='m',label='Potential Energy')
    axs[0].set_ylabel('Potential Energy')
    axs[0].set_xlabel('t')
    axs[0].set_ylim(-5*numOfParticles,5*numOfParticles)
    axs[0].grid()
    # Kinetic energy
    axs[1].plot(time, T,color='g',label='Kinetic Energy')
    axs[1].set_ylabel('Kinetic Energy')
    axs[1].set_xlabel('t')
    axs[1].set_ylim(0,5*numOfParticles)
    axs[1].grid()
    # Total energy
    axs[2].plot(time[0:-2], U+T ,color='b',label='Total Energy')
    axs[2].set_ylabel('Total Energy')
    axs[2].set_xlabel('t')
    axs[2].set_ylim(0,5*numOfParticles)
    axs[2].grid()
    # Set figure title
    if algorithm == "euler":
        fig.suptitle('Euler - {} particles, dt = {}'.format(numOfParticles,
                          timestep), size = 14)
    else:
        fig.suptitle('Velocity-Verlet - {} particles, dt = {}'.format(numOfParticles,
                          timestep), size = 14)
    plt.show()
    # Save figure
    if saveFigures == True:
        plt.savefig('{}.png'.format(name), dpi=300)

def pressure(MDS_dict):

    n = MDS_dict['num_t']
    n1 = MDS_dict['equilibrium_timestep']
    B = 119.8/MDS_dict['temp']
    rho = numOfParticles/(L**3)

    S = 0
    for ts in range(n1,n): # plus one to get last index
        for i, j in combinations(range(numOfParticles), 2):
            r,_,_,_ = getParticleDistance(i,j,ts)
            S = S + 24*(-2*(1/r)**11 + (1/r)**5) # S displays the sum over du/dr
    pp = 1 - (B/(6*numOfParticles)*(1/(n-n1)))*S
    P = pp*rho/B
    return P



################# Begin main program ########################

def main(MDS_dict):

    #######################
    ### Load Parameters ###
    #######################

    # Load parameters from dictionary into simulation
    algorithm = MDS_dict['algorithm']
    numOfTimesteps = MDS_dict['number_of_timesteps']
    timestep = MDS_dict['timestep']
    initParticles = MDS_dict['init_particles']
    latticeConstant = MDS_dict['lattice_constant']
    boxSize = MDS_dict['box_size']
    numOfParticles = MDS_dict['number_of_particles']
    numOfDimensions = MDS_dict['number_of_dimensions']
    bathTemperature = MDS_dict['bath_temperature']
    plotting = MDS_dict['plotting']
    plotCounter = MDS_dict['plot_counter']
    energyPlot = MDS_dict['energy_plot']
    saveFigures = MDS_dict['save_figures']

    pairCorrelation = MDS_dict['pair_correlation']
    # TODO load other observable flags

    #############################
    ### Simulation Parameters ###
    #############################

    # set simulation constants

    # number of timesteps stored before overwriting position, velocity, force data
    numOfStoredTimesteps = 100

    # number of unique particle pairs
    numberOfPairwiseInteractions = numOfParticles*(numOfParticles-1)/2

    # variable to track simulation equilibrium
    equilibriumTimestep = -1

    # number of timesteps before equilibration
    equilibrationTimer = 30

    ########################
    ### Simulation Setup ###
    ########################

    # Create n x d x 3 numpy array of floats "PC3T" to store n particles
    # P = particle, C = coordinate (x, y or z), 3 = (position,velocity, force), T = timestep
    # in d dimensions with 3 (position, velocity and force) parameters
    # The number of timesteps stored is specified by numOfStoredTimesteps.
    # Upon reaching this limit, the matrix will loop to the beginning and overwrite
    # old parameters.

    # For example, access particle 2's position in the y-direction at time_step=0 as
    # PC3T[1][1][0][0]
    # Access particle 1's momentum in the x-direction at time_step=1 as
    # PC3T[0][0][1][1]
    # Initialization later
    #PC3T = np.zeros((numOfParticles,numOfDimensions,3,numOfStoredTimesteps), dtype=float)

    # Initialize matrices to hold potential and kinetic energies, U and T respectively.
    # These arrays are of length numOfTimesteps as they track the entire simulation
    U = np.zeros((numOfTimesteps,1), dtype=float)
    T = np.zeros((numOfTimesteps,1), dtype=float)

    # Initialize data sturctures according to user-defined flags in main.py
    if pairCorrelation == True:
        # track particle distances to find pair correlation observable
        particleDistances = np.zeros((int(numberOfPairwiseInteractions),1), dtype=float)


    # set initial positions/velocities for all particles according to user selection
    if initParticles == 'fcc':
        numOfDimensions = 3 # Always use 3 dimensions (program hardcoded)
        # place particles on fcc lattice and initialize velocities with M-B distribution
        initialX, numOfParticles, boxSize = init_position(latticeConstant, numOfParticles, numOfDimensions)
        initialV = gaussVel(bathTemperature, numOfParticles, numOfDimensions)
        PC3T = np.zeros((numOfParticles,numOfDimensions,3,numOfStoredTimesteps), dtype=float)
        print(initialX.shape)
        print(PC3T.shape)
        PC3T[:][:][0][0] = initialX
        PC3T[:][:][1][0] = initialV

    elif initParticles == 'random':
        raise Exception('Not coded yet')

    # TODO modify debug modes to work with cleaned-up code. We cannot use old set/get functions
    elif initParticles == 'debug':
        raise Exception('Not coded yet')

    if plotting == True:
        # set up scatter plot for displaying particle motion
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatterPoints = [] # hold scatter plot data for each particle
        colours = ['b','g','r','c','m','y','k','w']
        for p in range(0,numOfParticles):
            colour = colours[p%7]
            scatterPoints.append(ax.scatter(PC3T[p,0,0,0],PC3T[p,1,0,0],PC3T[p,2,0,0],color=colour))
        ax.set_xlim((0,boxSize))
        ax.set_ylim((0,boxSize))
        ax.set_zlim((0,boxSize))
        plt.ion()
        plt.show()
        plt.pause(0.01)

    # TODO figure out a solution for the wrap around of timestep index

    # use iteration in range(num_t) and iterateCoordinates for euler method
    # use iteration in range(num_t-1) due to need of "future" force in Verlet algorithm

    ##################
    ### Simulation ###
    ##################

    if algorithm == "verlet":
        # need to initialize force for Verlet algorithm
        particleDistances = getParticleDistanceArray(PC3T[:,:,0,0],numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,2,0], U[0] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)

    for j in range(numOfTimesteps):
        i = j%(numOfStoredTimesteps) # overwrite parameter matrix

        # First calculate kinetic energy at the given timestep
        T[j] = getKineticEnergy(PC3T[:,:,1,i],numOfParticles,numOfDimensions)

        # TODO investigate how to best establish when equilibrium has been reached
        # Equilibrate according to simulation constants set above
        if i%equilibrationTimer == equilibrationTimer-1:
            # check if at equilibrium, note that 119.8 converts natural energy scale to Kelvin (see week 1 notes)
            # see if temperature agrees to desired temperature within <1> Kelvin. If it does, don't equilibrate
            if equilibriumTimestep == -1:
                if abs(float(T[j]) / (numOfParticles-1) / (3/2) * 119.8 - bathTemperature) > 1:
                    # we need to equilibrate
                    # scale all particle velocities
                    scalingConstant = getVelocityScalingFactor(T[j],bathTemperature)
                    PC3T[:,:,1,i] = np.multiply(PC3T[:,:,1,i],scalingConstant)

                    # recalculate present kinetic energy
                    print('Rescaling from temperature: ',float(T[j]) / (numOfParticles-1) / (3/2) * 119.8)
                    T[j] = getKineticEnergy(PC3T[:,:,1,i],numOfParticles,numOfDimensions)
                else:
                    # we can use equilibrium_timestep as the beginning index for calculating observables
                    equilibriumTimestep = j
                    MDS_dict['equilibrium_timestep'] = equilibriumTimestep
                    print('Rescaled temperature: ',float(T[j]) / (numOfParticles-1) / (3/2) * 119.8)

                    if plotting == True:
                        plt.title('Equilibrium reached')

                    # zero observables after equilibrium os reached
                    # TODO only start tracking observables after equilibrium
                    #particleDistances = np.multiply(particleDistances,0)

        ### Measure Observables ###
        # track observables here; this can also be done using variables
        # calculated in the time evolution block to avoid double-calculations

        ### Time Evolution ###
        if algorithm == "euler":
            # Euler time-evolution

            # Calculate inter-particle distances at current timestep
            particleDistances = getParticleDistanceArray(PC3T[:,:,0,i],numOfParticles,numOfDimensions,boxSize)
            # Calculate forces at current timestep
            PC3T[:,:,2,i],U[j] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)
            if i+1 < numOfStoredTimesteps:
                PC3T[:,:,0,i+1] = iterateCoordinates_Euler(PC3T[:,:,0,i],PC3T[:,:,1,i],numOfParticles,numOfDimensions,timestep)
                PC3T[:,:,1,i+1] = iterateVelocities_Euler(PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep)
            else:
                PC3T[:,:,0,0] = iterateCoordinates_Euler(PC3T[:,:,0,i],PC3T[:,:,1,i],numOfParticles,numOfDimensions,timestep)
                PC3T[:,:,1,0] = iterateVelocities_Euler(PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep)
        else:
            # Verlet time-evolution
            # Force at time = 0 is calculated before evolving for Verlet method
            if i+1 < numOfStoredTimesteps:
                PC3T[:,:,0,i+1] = iterateCoordinates_Verlet(PC3T[:,:,0,i],PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep)
                # Must evaluate force at next timestep to iterate velocities
                # Calculate inter-particle distances at current timestep
                particleDistances = getParticleDistanceArray(PC3T[:,:,0,i+1],numOfParticles,numOfDimensions,boxSize)
                # Calculate forces at current timestep
                PC3T[:,:,2,i+1],U[j+1] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)
                PC3T[:,:,1,i+1] = iterateVelocities_Verlet(PC3T[:,:,1,i],PC3T[:,:,2,i],PC3T[:,:,2,i+1],numOfParticles,numOfDimensions,timestep)
            else:
                PC3T[:,:,0,0] = iterateCoordinates_Verlet(PC3T[:,:,0,i],PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep)
                # Evaluate force at next timestep, as above
                particleDistances = getParticleDistanceArray(PC3T[:,:,0,0],numOfParticles,numOfDimensions,boxSize)
                PC3T[:,:,2,0],U[j+1] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)

                PC3T[:,:,1,0] = iterateVelocities_Verlet(PC3T[:,:,1,i],PC3T[:,:,2,i],PC3T[:,:,2,0],numOfParticles,numOfDimensions,timestep)

        if plotting == True:
            if j%plotCounter == 0:
                for p in range(len(scatterPoints)):
                    scatterPoints[p].remove()
                    colour = colours[p%7]
                    scatterPoints.append(ax.scatter(PC3T[p,0,0,i],PC3T[p,1,0,i],PC3T[p,2,0,i],color=colour))
                plt.pause(0.000005)

    ########################
    ### Post-proecessing ###
    ########################

    if energyPlot == True:
        plotEnergy(algorithm,numOfParticles,numOfTimesteps,timestep,saveFigures,U,T)

    ### Compute Observables ###

    '''
    P = pressure(MDS_dict)
    print('Pressure of the System: P = {}'.format(P))

    numOfBins = 100 # number of bins for pair correlation histogram

    # create histogram for pair correlation function
    # the histogram of the data
    plt.figure(2)
    particleDistances = particleDistances/(num_t-1-equilibrium_timestep)
    plt.ion()
    n, bins, patches = plt.hist(particleDistances, numOfBins, facecolor='g')

    # Calculate pair correlation function for each bin
    pairCorrelation = np.zeros((numOfBins,1),dtype=float)
    for i in range(numOfBins):
        pairCorrelation[i] = 2*L**3/numOfParticles/(numOfParticles-1)*n[i]/4/np.pi/(bins[i+1])**2/(bins[i+1]-bins[i])

    plt.figure(3)
    plt.plot(bins[0:-1],pairCorrelation)
    plt.ylabel('g(r)')
    plt.xlabel(r'$r/\sigma$')
    plt.show()
    plt.pause(30)
    '''



    #print(sum(particleDistances[equilibrium_timestep:-2])/(num_t-2-equilibrium_timestep))


    return MDS_dict
