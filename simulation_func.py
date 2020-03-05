from mpl_toolkits.mplot3d import Axes3D # important for 3d scatter plot
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
from scipy.optimize import fsolve
import os
=======
from math import sqrt
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e

##### Parameters #####

# Lennard-Jones parameters Argon
#eps = 119.8 # k_b
#sigma = 3.405 # Angstrom
#mass = 39.948*1.66e-27 #Kg


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

def getKineticEnergy(ts):
    KE = 0
    for i in range(numOfParticles):
        KE = KE + 0.5*(getPXvel(i,ts)**2 + getPYvel(i,ts)**2 + getPZvel(i,ts)**2)
    T[ts] = KE


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

# Euler algorithm
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

# Verlet algorithm
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

# Verlet algorithm
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

<<<<<<< HEAD
def init_position(ts):
    # Initializes particles into fcc lattice
    # for a given n = L/a, an fcc compatible number of particles is
    # N = (n + 1)**3 + 3*(n+1)*n**2

    global PC3T
    global numOfParticles
    global L

    # Find a compatible n
    func = lambda x : numOfParticles - (x + 1)**3 - 3*x**3
    n = int(np.round(fsolve(func, 1.1)))
    # Compute L and exact N
    L = a*n + a # +a to avoid putting particles on boundary
    numOfParticles = (n + 1)**3 + 3*(n+1)*n**2

    # Reinitialize PC3T
    PC3T = np.zeros((numOfParticles,numOfDimensions,3,num_t), dtype=float)
    
    # Put particles on cubic lattice
    for i in range(n+1): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n+1): # loop over x
                p = i*(n+1)**2 + j*(n+1) + k #particle index, (n+1)**3 total on the cubic lattice
                setPXcoord(k*a + a/2,p,ts) #a/2 avoids putting particles on boundary
                setPYcoord(j*a + a/2,p,ts)
                setPZcoord(i*a + a/2,p,ts)
                #print("x:{}\n y:{}\n z:{}".format(k*a + a/2, j*a + a/2, i*a + a/2))
    # Put particles on every xz face
    for i in range(n): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n): # loop over x
                p = (n+1)**3 + i*n*(n+1) + j*n + k #particle index, (n+1)*n**2 on x faces
                setPXcoord(k*a + a,p,ts)
                setPYcoord(j*a + a/2,p,ts)
                setPZcoord(i*a + a,p,ts)
    # Put particles on every yz face
    for i in range(n): # loop over z
        for j in range(n): # loop over y
            for k in range(n+1): # loop over x
                p = (n+1)**3 + (n+1)*n**2 + i*n*(n+1) + j*(n+1) + k #particle index, (n+1)*n**2 on yz faces
                setPXcoord(k*a + a/2,p,ts)
                setPYcoord(j*a + a,p,ts)
                setPZcoord(i*a + a,p,ts)
    # Put particles on every xy face
    for i in range(n+1): # loop over z
        for j in range(n): # loop over y
            for k in range(n): # loop over x
                p = (n+1)**3 + 2*(n+1)*n**2 + i*n*n + j*n + k #particle index, (n+1)*n**2 on xy faces
                setPXcoord(k*a + a,p,ts)
                setPYcoord(j*a + a,p,ts)
                setPZcoord(i*a + a/2,p,ts)
=======
def gaussVel(T,N):
    """
    Function to initialize particles velocity components according to Gaussian distribution.
    A function called velocityAVG is called for each velocity component array to ensure the avg velocity to be zero.
    
    Parameters
    -------------------------
    T   Temperture in Kelvin
    N   Number of atoms

    Return
    -------------------------
    vx,vx,vz as numpy arrarys of length N holding velocity components for N atoms
    the average of each velocity component array is zero.
    
    """ 
    
    mu, sigma = 0, sqrt((119.8*T)/2) # mean is 0 and standard deviation in Kelvin
    vx = np.random.normal(mu, sigma, N)
    vy = np.random.normal(mu, sigma, N)
    vz = np.random.normal(mu, sigma, N)
    
    vx = velocityAVG(vx)
    vy = velocityAVG(vy)
    vz = velocityAVG(vz)

    return vx,vy,vz

def velocityAVG(a):
    '''
    Function to set the mean value of the component velocities close to zero by adjusting each velocity
    
    Parameters
    -------------------------
    a   np.array holding velocity components 
    
    Return
    -------------------------
    vx,vx,vz as numpy arrarys of length N holding velocity components for N atoms
    the average of each velocity component array is zero.      
    '''
    x = np.mean(a)
    while not -0.01 < x < 0.01:
        a = a -(np.sum(a)/len(a))
        x = np.mean(a)
    return a
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e

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

    # this variable is now included in main.py
    maxInitialSpeedParameter = 1000

    dim_components = np.zeros((numOfDimensions,1),dtype=float)

    for i in range(0,numOfParticles):
        # first generate positions in each dimension
        dim_components = np.random.rand(numOfDimensions,1)*L # scale to sizeOfBox
        dim_components[0] = dim_components[0]/numOfParticles/1.1 + i/L # slice in dimension 1 to separate particles initially (1.1 gives space between "sliced" regions)
        # pass vector to method to fill parameter matrix:
        # arguments: random vector, particle number, coord/vel, timestep
        addToParameterMatrix(dim_components,i,0,ts)

        # next generate velocities in each dimension, limited according to maxInitialSpeedParameter above
        dim_components = np.random.rand(numOfDimensions,1)/np.sqrt(numOfDimensions)*L/timestep/maxInitialSpeedParameter
        # scale velocities to be either positive or negative
        dim_components = dim_components*2 - L/timestep/maxInitialSpeedParameter/np.sqrt(numOfDimensions)
        addToParameterMatrix(dim_components,i,1,ts)

def addToParameterMatrix(dim_components,pnum,xv,ts):
    # treat this  as Protected method
    #
    # function called by initializePartciels() to load randomly generated initial
    # positions/velocities into parameter matrix

    if xv == 0:
        # load positions
        for d in range(0,len(dim_components)):
            setPncoord(dim_components[d],d,pnum,ts)
    elif xv == 1:
        # load velocities
        for d in range(0,len(dim_components)):
            setPnvel(dim_components[d],d,pnum,ts)

def plotEnergy(MDS_dict):

    ##### Plots #########
    # Euler
    if MDS_dict['euler'] == True:
        name = 'MDS_{}p_{}'.format(MDS_dict['numOfParticles'],'euler')
        fig = plt.figure(name)
    # Verlet
    else:
        name = 'MDS_{}p_{}'.format(MDS_dict['numOfParticles'],'verlet')
        fig = plt.figure(name)

    plt.ioff()
    time = np.arange(0, num_t*timestep, timestep)
    a = fig.subplots(2,2)
    a[0][0].plot(time[0:-2],particleDistances[0:-2],color='r',label = 'Inter-particle distance')
    a[0][0].set_ylabel('Particle Distance')
    a[0][0].set_xlabel('time')
    a[0][1].plot(time[0:-2], U[0:-2],color='m',label='Potential Energy')
    a[0][1].set_ylabel('Potential Energy')
    a[0][1].set_xlabel('time')
    a[1][0].plot(time[0:-2], T[0:-2],color='g',label='Kinetic Energy')
    a[1][0].set_ylabel('Kinetic Energy')
    a[1][0].set_xlabel('time')
    a[1][1].plot(time[0:-2], (U[0:-2]+T[0:-2])/numOfParticles,color='b',label='Total Energy')
    a[1][1].set_ylabel('Total Energy')
    a[1][1].set_xlabel('time')
    #a[1][1].set_ylim(0,1000)
    a[0][0].grid()
    a[0][1].grid()
    a[1][0].grid()
    a[1][1].grid()
    # Euler
    if MDS_dict['euler'] == True:
        fig.suptitle('Euler - {} particles, dt = {}'.format(MDS_dict['numOfParticles'],
                          MDS_dict['timestep']), size = 14)
    # Verlet
    else:
        fig.suptitle('Velocity-Verlet - {} particles, dt = {}'.format(MDS_dict['numOfParticles'],
                          MDS_dict['timestep']), size = 14)
    plt.show()

    if MDS_dict['save_fig'] == True:
#        path = os.getcwd()
        plt.savefig('{}.png'.format(name), dpi=150)

    #print('Initial Energy:', E[0],'\nFinal Energy', E[-2])

<<<<<<< HEAD
def dictTester(D):

=======
def printer():

    print('print from function: {}'.format(numOfParticles))

def dictTester(D):
    """Func to check the input parameters saved in the dictionary MDS_dict
    
    Parameters
    -------------------------
    D   MDS_dict

    Raise
    -------------------------
    ValueError: If verlet and euler method are both chosen to be True or False
    ValueError: If the debug_2 mode is chosen with particle number differing from 2
    
    """
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e
#    try:
#        D['euler'] or D['verlet'] == False
#    except ValueError as err:
##        sys.stderr.write("{}: Cannot open file {}: {}\ n".format(sys.argv[0], filename, err ))
#        print('{} Only one method at a time coputable.'.format(err))
#        sys.exit(1)
    if (D['euler'] and D['verlet']) or (not D['euler'] and not D['verlet']):
      raise ValueError('Only one method at a time computable.')

    if D['init_particles'] == 'debug_2' and D['numOfParticles'] != 2 :
      raise ValueError('Debug_2 assigns positions for 2 particles. \
                       Use another initialisation method for more particles.')



#
#    MDS_parameters = {
#'euler' :           False,
#'verlet' :          True,
#'boxSize_L' :       5,
#'numOfParticles' :  2,
#'numOfDimensions' : 3,
#'num_t' :           200,
#'timestep' :        0.001,
#'plotting':         True,
#'plot_counter' :    1,
#'energyPlot' :      True,
#'init_particles' :  'debug_2'
#}

################# Begin main program ########################

def main(MDS_dict):

    dictTester(MDS_dict)

    # setting of the used parameters stored in dictionary
    global L
    global a
    L = MDS_dict['boxSize_L'] # size of each periodic unit cell L (in units of sigma)
    a = MDS_dict['latticeConst'] # lattice constant
    # init_position will overwrite L when called

    global numOfParticles
    global numOfDimensions
    numOfParticles = MDS_dict['numOfParticles'] # 2 particles
    numOfDimensions = MDS_dict['numOfDimensions'] # 3D
<<<<<<< HEAD

=======
    
    printer()
    print('print from main: {}'.format(numOfParticles))
    
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e
    #Total_time = 500*0.001
    global num_t
    global timestep
    num_t = MDS_dict['num_t']
    timestep = MDS_dict['timestep']

    global plotting
    global plot_counter
    plotting = MDS_dict['plotting']
    plot_counter = MDS_dict['plot_counter']

    init_particles = MDS_dict['init_particles']


    # Create n x d x 3 numpy array of floats "PC3T" to store n particles
    # P = particle, C = coordinate (x, y or z), 3 = (position,velocity, force), T = timestep
    # in d dimensions with 3 (position, velocity and force) parameters, and num_t timesteps stored.

    # For example, access particle 2's position in the y-direction at time_step=0 as
    # PC3T[1][1][0][0]
    # Access particle 1's momentum in the x-direction at time_step=1 as
    # PC3T[0][0][1][1]

    global PC3T
    PC3T = np.zeros((numOfParticles,numOfDimensions,3,num_t), dtype=float)

    # Initialize potential energy matrix U and total energy E
    # for tracking potential and total energy at each timestep
    global U
    global E
    global T
    U = np.zeros((num_t,1), dtype=float)
    E = np.zeros((num_t,1), dtype=float)
    T = np.zeros((num_t,1), dtype=float)

    ##### Set initial positions/velocities for all particles ####
    if init_particles == 'fcc':
        init_position(0)
    elif init_particles == 'random':
        initializeParticles(0)
    elif init_particles == 'debug_2':
        # Particle 1
        setPXcoord(L/2-0.5,0,0)
        setPYcoord(L/2,0,0)
        setPZcoord(0,0,0)
        setPXvel(0.5,0,0)
        setPYvel(0,0,0)
        setPZvel(0,0,0)
        # Particle 2
        setPXcoord(L/2+0.5,1,0)
        setPYcoord(L/2,1,0)
        setPZcoord(0,1,0)
        setPXvel(-0.5,1,0)
        setPYvel(0,1,0)
        setPZvel(0,1,0)
        # Particle 3
#        setPXcoord(3*L/4,2,0)
#        setPYcoord(3*L/4,2,0)
#        setPZcoord(L/2,0,0)
#        setPXvel(0,1,0)
#        setPYvel(0,1,0)
#        setPZvel(0,1,0)
    else:
        print('choose appropriate string for particle initialisation: \n\
<<<<<<< HEAD
              {} \n{}\n{}'.format('fcc','random','debug_3'))


=======
              {} \n{}\n{}'.format('fcc','random','debug_2'))
        
        
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e
    ##### Simulation #####

    if plotting == True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatterPoints = [] # hold scatter plot data for each particle
        colours = ['b','g','r','c','m','y','k','w']
        for p in range(0,numOfParticles):
            colour = colours[p%7]
            scatterPoints.append(ax.scatter(getPXcoord(p,0),getPYcoord(p,0),getPZcoord(p,0),color=colour))
        ax.set_xlim((0,L))
        ax.set_ylim((0,L))
        ax.set_zlim((0,L))
        plt.ion()
        plt.show()
        plt.pause(0.01)

    # vp1 and vp2 to keep track of the distances
    vp1 = np.zeros(num_t)
    vp2 = np.zeros(num_t)

    # keep track of particle distance, potential energy, and kinetic energy
    global particleDistances
    particleDistances = np.zeros(num_t)

    # TODO figure out a solution for the wrap around of timestep index

    # use iteration in range(num_t) and iterateCoordinates for euler method
    # use iteration in range(num_t-1) due to need of "future" force in Verlet algorithm

    # TODO find a better implementation of this

    # REMOVE THIS FOR EULER METHOD
    # Get forces for initial positions
    if MDS_dict['verlet'] == True:
        getForce(0)

    for j in range(num_t-1):
        i = j%(num_t-1) # don't go over indices of PC3

        #Euler
        if MDS_dict['euler'] == True:
            iterateCoordinates(i)
            iterateVelocities(i)

        #Velocity-Verlet
        if MDS_dict['verlet'] == True:
            iterateCoordinates_Verlet(i)
            iterateVelocities_Verlet(i)

        #getTotalEnergy(j)
        getKineticEnergy(j)

        if plotting == True:
            if j%plot_counter == 0:
                for p in range(len(scatterPoints)):
                    scatterPoints[p].remove()
                    colour = colours[p%7]
                    scatterPoints[p] = ax.scatter(getPXcoord(p,i),getPYcoord(p,i),getPZcoord(p,i),color=colour)
                plt.pause(0.000005)

#        p1.remove()
#        p2.remove()
#        p1 = ax.scatter(getPXcoord(0,i),getPYcoord(0,i), getPZcoord(0,i),color='r')
#        p2 = ax.scatter(getPXcoord(1,i),getPYcoord(1,i), getPZcoord(0,i),color='b')

        vp1[i] = getPXvel(0,i)
        vp2[i] = getPXvel(1,i)

        particleDistances[i],a,b,c = getParticleDistance(0,1,i)
        #time.sleep(0.05)

    if MDS_dict['energyPlot'] == True:
        plotEnergy(MDS_dict)

    return MDS_dict
