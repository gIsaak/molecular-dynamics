from mpl_toolkits.mplot3d import Axes3D # important for 3d scatter plot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import os
from itertools import combinations
##### Parameters #####

# Lennard-Jones parameters Argon
#eps = 119.8 # k_b
#sigma = 3.405 # Angstrom
#mass = 39.948*1.66e-27 #Kg

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
    particleDistances = np.zeros((int((numOfParticles**2-numOfParticles)/2),1+numOfDimensions),dtype=float)
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
    return newV

def init_position(a,numOfParticles,numOfDimensions,density):
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
    func = lambda x : numOfParticles - (x + 1)**3 - 3*x**3 # Find a compatible n
    n = int(np.round(fsolve(func, 1.1))) #1.1 initial guess

    # Compute L and exact N
    L = a*n + a # +a to avoid putting particles on boundary
    newNumOfParticles = (n + 1)**3 + 3*(n+1)*n**2

    # if statement to adjust lattice constant in case density is specified
    if density != 0:
        a = ((newNumOfParticles/density)**(1/3))/(n+1)
        L = a*n + a # +a to avoid putting particles on boundary
    print('Density: {}\nNumber N: {}\nLattice const: {}'.format(newNumOfParticles/(L**3),newNumOfParticles,a))

    newX = np.zeros((newNumOfParticles,numOfDimensions), dtype=float)
    p = 0 # particle index
    # Put particles on cubic lattice
    for i in range(n+1): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n+1): # loop over xe
                newX[p,0] = k*a + a/2  #a/2 avoids putting particles on boundary
                newX[p,1] = j*a + a/2
                newX[p,2] = i*a + a/2
                p += 1
    # Put particles on xz faces
    for i in range(n): # loop over z
        for j in range(n+1): # loop over y
            for k in range(n): # loop over x
                newX[p,0] = k*a + a
                newX[p,1] = j*a + a/2
                newX[p,2] = i*a + a
                p += 1
    # Put particles on yz faces
    for i in range(n): # loop over z
        for j in range(n): # loop over y
            for k in range(n+1): # loop over x
                newX[p,0] = k*a + a/2
                newX[p,1] = j*a + a
                newX[p,2] = i*a + a
                p += 1
    # Put particles on xy faces
    for i in range(n+1): # loop over z
        for j in range(n): # loop over y
            for k in range(n): # loop over x
                newX[p,0] = k*a + a
                newX[p,1] = j*a + a
                newX[p,2] = i*a + a/2
                p += 1
    return newX, newNumOfParticles, L

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
    return newV


def getVelocityScalingFactor(KE,bathTemperature,numOfParticles):
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
def plotEnergy(numOfParticles,numOfTimesteps,timestep,saveFigures,U,T):
    '''
    Creates plot for Kinetic Energy, Potential energy and Total energy.
    Plot title/name/saving/etc is handled given the parameters set in main.py

    The plots will show all energies at all steps, including before the
    desired equilibrium temperature has been reached.
    '''
    # Figure name
    name = 'Eplot_N_{}'.format(numOfParticles)
    fig = plt.figure(name)
    plt.ioff()
    time = np.arange(0, numOfTimesteps*timestep, timestep)
    plt.plot(time, U[0:-1],color='m',label='Potential Energy')
    plt.plot(time, T[0:-1],color='g',label='Kinetic Energy')
    plt.plot(time, U[0:-1]+T[0:-1] ,color='b',label='Total Energy')
    plt.xlabel('t')
    plt.ylabel('Energy')
    plt.legend()
    plt.grid()
    fig.suptitle('{} particles, dt = {}'.format(numOfParticles,timestep), size = 14)
    plt.show()
    # Save figure
    if saveFigures == True:
        plt.savefig('{}.png'.format(name), dpi=300)

def autocorr(data, plot=True):
    '''
    Computes autocorrelation function for an observable

    Autocorrelation function as calculated in the lecutre notes (week 5)
    The cutoff for t is given by t = sqrt(N_sim)
    Autocorrelation length tau is computed fitting the curve to
    a*exp(-t/tau)

    Parameters
    ----------

    data: numpy array
        value of observable at each timestep

    Returns:
    --------
    chi: numpy array
        autocorrelation function
    tau: float
        autocorrentaion length
    '''
    N_sim = int(data.size) #number of simulation steps
    chi_length = int(round(np.sqrt(N_sim))) #cutoff at sqrt(N_sim)
    chi = np.zeros(shape=(chi_length,))
    for t in range(chi_length):
        a, b, c = 0, 0, 0
        for n in range(N_sim - t):
            a = a + data[n]*data[n+t]
            b = b + data[n]
            c = c + data[n+t]
        chi[t] = a/(N_sim - t) - b*c/(N_sim - t)**2
    # Fit
    xdata = np.arange(chi_length)
    def func(x, tau, a):
        return np.exp(-x/tau)*a
    param, _ = curve_fit(func, xdata, chi)
    tau = param[0]
    if plot == True:
        plt.plot(xdata, chi, 'b-', label='data')
        # Fit
        plt.plot(xdata, func(xdata, *param), 'r--', label=r'Autocorrelation $\tau$=%5.3f, a=%5.3f' % tuple(param))
        plt.xlabel('t')
        plt.ylabel(r'$\chi$')
        plt.title(r'N = %d' % N_sim)
        plt.legend()
        plt.show()
    return chi, tau

def averageAndError(obs,autoCorr=False, plotAutoCorr=False):
    '''
    Function to computer and observable average and error
    - autoCorr = False (default) assumes statistically independent data
    - autoCorr = True            rescales with correlation time tau

    Accepts: obs: array of observable values at each time-step
             numOfTimesteps: number of observable values
    Returns: avg: avergae of observable
             sigma: observable error
    '''
    obs2 = np.square(obs)
    avg = np.mean(obs)
    avg2 = np.mean(obs2)
    N = obs.size
    sigma = np.sqrt((avg2 - avg**2)/(N - 1))
    if autoCorr == True:
        _, tau = autocorr(obs, plotAutoCorr)
        sigma = sigma*np.sqrt(2*tau*(N - 1)/N)
    return avg, sigma

def pcf_count(distances, numOfBins, boxSize):
    '''
    Function to create histcount to compute PCF (Pair Correlation  Function)
    Accepts: distances: array of particle-particle distances
             numOfBins: number of bins in hsitogram
             boxSize: simultion box size
    Returns: histCount array
    '''
    histCount = np.zeros(shape=(numOfBins,))
    dr = boxSize/numOfBins
    for j in range(distances.size):
        i = int(distances[j]//dr)
        histCount[i] += 1
    return histCount

def pcf_plot(histCount, numOfParticles, numOfTimesteps, numOfBins, boxSize, saveFigures = False):
    '''
    Function to compute and plot PCF (Pair Correlation Function)
    Accepts: histCount: array pcf histcounts at each timestep shape=(numOfBins,numOfTimesteps)
             numOfParticles: simulation particles
             numOfTimesteps: number of timesteps used to compute PCF
             numOfBins: number of bins in hsitogram
             boxSize: simultion box size
    Returns: pcf array
    '''
    pcf = np.zeros(shape=(numOfBins,), dtype= float)
    pcfErr = np.zeros(shape=(numOfBins,), dtype= float)
    r = np.zeros(shape=(numOfBins,), dtype= float)

    numOfParticlePairs = numOfParticles*(numOfParticles-1)/2
    V = boxSize**3
    dr = boxSize/numOfBins
    for j in range(pcf.size):
        r[j] = (j + 0.5)*dr
        rawPCF, rawErrPCF = averageAndError(histCount[:,j])
        pcf[j]= rawPCF *V/numOfParticlePairs/(4*np.pi*dr*r[j]**2)
        pcfErr[j] = rawErrPCF *V/numOfParticlePairs/(4*np.pi*dr*r[j]**2)
    # Plot
    name = 'pcfN_{}'.format(numOfParticles)
    fig = plt.figure(name)
    plt.ioff()
    plt.errorbar(r, pcf, pcfErr, label='Pair Correlation Function')
    plt.xlabel('r')
    plt.ylabel('PCF')
    plt.grid()
    fig.suptitle('Pair correlation function', size = 14)
    plt.show()
    # Save figure
    if saveFigures == True:
        plt.savefig('{}.png'.format(name), dpi=300)
    return pcf, pcfErr


def pressure(particleDistances):
    '''
    Function to compute pressure at time t
    Accepts: particleDistances: array of particle distances
             numOfParticles: number of particles in the simulation
             bathTemperature: temperature
    Returns: P: pressure
    '''
    P = 0
    S = 0
    for i in range(particleDistances[:,0].size):
        r = particleDistances[i,0]
        invr6 = (1/r)**6 #precomputes (1/r)**6
        S = S + 12* (-2*invr6**2  + invr6)
    P = S
    return P

def diffusion(numOfParticles,numOfDimensions,initialX,PC3T):
    '''
    Function to gather position data of each particle to calculate diffusion
    Accepts: numOfParticles: number of particles in the simulation
             numOfDimensions: number of dimensions
             initialX: initial positions of particles
             PC3T[:,:,0,i]: all particles x,y,z positions at timestep i
    Returns: diff: array with mean square distance moved by all particles
             
    '''
    c_diff = np.zeros(shape=(numOfParticles,2))
    diff = np.zeros(shape=(numOfParticles))
    #for g in range(numOfParticles):
    for d in range(numOfDimensions):
        c_diff[:,0] = c_diff[:,0] + (PC3T[:,d])**2
        c_diff[:,1] = c_diff[:,1] + (initialX[:,d])**2
    diff = (np.sqrt(c_diff[:,0])-np.sqrt(c_diff[:,1]))**2
    return np.mean(diff)
            
################# Begin main program ########################

def main(MDS_dict):

    #######################
    ### Load Parameters ###
    #######################

    # Load parameters from dictionary into simulation
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
    density = MDS_dict['density']

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
    equilibrationTimer = 20

    # initialize memory for plotting
    scatterPoints = []
    colours = ['b','g','r','c','m','y','k','w']

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
    U = np.zeros((numOfTimesteps+1,1), dtype=float)
    T = np.zeros((numOfTimesteps+1,1), dtype=float)
    particleDistances = np.zeros((int(numberOfPairwiseInteractions),1), dtype=float)


    # set initial positions/velocities for all particles according to user selection
    if initParticles == 'fcc':
        numOfDimensions = 3 # Always use 3 dimensions (program hardcoded)
        # place particles on fcc lattice and initialize velocities with M-B distribution
        initialX, numOfParticles, boxSize = init_position(latticeConstant, numOfParticles, numOfDimensions,density)
        initialV = gaussVel(bathTemperature, numOfParticles, numOfDimensions)
        PC3T = np.zeros((numOfParticles,numOfDimensions,3,numOfStoredTimesteps), dtype=float)
        PC3T[:,:,0,0] = initialX
        PC3T[:,:,1,0] = initialV
        particleDistances = getParticleDistanceArray(PC3T[:,:,0,0],numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,2,0], U[0] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)
    elif initParticles == 'random':
        raise Exception('Not coded yet')
    elif initParticles == 'debug':
        raise Exception('Not coded yet')

    if plotting == True:
        # set up scatter plot for displaying particle motion
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for p in range(numOfParticles):
            colour = colours[p%7]
            scatterPoints.append(ax.scatter(PC3T[p,0,0,0],PC3T[p,1,0,0],PC3T[p,2,0,0],color=colour))
        ax.set_xlim((0,boxSize))
        ax.set_ylim((0,boxSize))
        ax.set_zlim((0,boxSize))
        plt.ion()
        plt.show()
        plt.pause(0.01)


    ##################
    ### Simulation ###
    ##################

    ### EQUILIBRATION STEPS ###
    equilibrium = False
    j = 0
    T_equilibration = 0
    while equilibrium == False:
        i = j%(numOfStoredTimesteps) # overwrite parameter matrix
        i_next = (j+1)%(numOfStoredTimesteps)
        # Evolution
        PC3T[:,:,0,i_next] = iterateCoordinates_Verlet(PC3T[:,:,0,i],PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep,boxSize)
        particleDistances = getParticleDistanceArray(PC3T[:,:,0,i_next],numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,2,i_next], _ = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,1,i_next] = iterateVelocities_Verlet(PC3T[:,:,1,i],PC3T[:,:,2,i],PC3T[:,:,2,i_next],numOfParticles,numOfDimensions,timestep)
        # scatter plot
        if plotting == True:
            if j%plotCounter == 0:
                for p in range(len(scatterPoints)):
                    scatterPoints.pop(0).remove()
                    colour = colours[p%7]
                    scatterPoints.append(ax.scatter(PC3T[p,0,0,i],PC3T[p,1,0,i],PC3T[p,2,0,i],color=colour))
                plt.pause(0.000005)
        # Equilibrium control
        T_equilibration = getKineticEnergy(PC3T[:,:,1,i],numOfParticles,numOfDimensions)
        # check if at equilibrium every equilibrationTimer steps
        # note that 119.8 converts natural energy scale to Kelvin (see week 1 notes)
        # see if temperature agrees to desired temperature within <1> Kelvin. If it does, don't equilibrate
        if i%equilibrationTimer == equilibrationTimer-1:
            if abs(float(T_equilibration) / (numOfParticles-1) / (3/2) * 119.8 - bathTemperature) > 1:
                scalingConstant = getVelocityScalingFactor(T_equilibration,bathTemperature,numOfParticles)
                print('Scaling constant:', scalingConstant)
                PC3T[:,:,1,i_next] = np.multiply(PC3T[:,:,1,i_next],scalingConstant) #rescale
                T_equilibration = getKineticEnergy(PC3T[:,:,1,i_next],numOfParticles,numOfDimensions) #update kinetic energy
                print('Rescaling from temperature: ',float(T_equilibration) / (numOfParticles-1) / (3/2) * 119.8)
            else:
                equilibrium = True
                # Set this configuration as initial configuration
                PC3T[:,:,:,0] = PC3T[:,:,:,i_next]
                print('Rescaled temperature: ',float(T_equilibration) / (numOfParticles-1) / (3/2) * 119.8)
                print('From initial temperature: {}'.format(bathTemperature))
                if plotting == True:
                    plt.title('Equilibrium reached')
        j += 1
    print('{} steps to equilibrium'.format(j))



    # Initialize observables
    numOfBins = int(round(np.sqrt(numOfTimesteps)))
    pcfCount = np.zeros(shape=(numOfTimesteps,numOfBins))
    P = np.zeros(shape=(numOfTimesteps,))
    D = np.zeros(shape=(numOfTimesteps))
    D_avg = np.zeros(shape=(numOfTimesteps))
    ### BEGINNING OF BIG LOOP ###
    for j in range(numOfTimesteps):
        i = j%(numOfStoredTimesteps) # overwrite parameter matrix
        i_next = (j+1)%(numOfStoredTimesteps)

        # First calculate kinetic energy at the given timestep
        T[j] = getKineticEnergy(PC3T[:,:,1,i],numOfParticles,numOfDimensions)
        #Evolution
        PC3T[:,:,0,i_next] = iterateCoordinates_Verlet(PC3T[:,:,0,i],PC3T[:,:,1,i],PC3T[:,:,2,i],numOfParticles,numOfDimensions,timestep,boxSize)
        particleDistances = getParticleDistanceArray(PC3T[:,:,0,i_next],numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,2,i_next], U[j+1] = getForceAndPotentialEnergy(particleDistances,numOfParticles,numOfDimensions,boxSize)
        PC3T[:,:,1,i_next] = iterateVelocities_Verlet(PC3T[:,:,1,i],PC3T[:,:,2,i],PC3T[:,:,2,i_next],numOfParticles,numOfDimensions,timestep)
        if plotting == True:
            if j%plotCounter == 0:
                for p in range(len(scatterPoints)):
                    scatterPoints.pop(0).remove()
                    colour = colours[p%7]
                    scatterPoints.append(ax.scatter(PC3T[p,0,0,i],PC3T[p,1,0,i],PC3T[p,2,0,i],color=colour))
                plt.pause(0.000005)
        #Observables
        pcfCount[j,:] = pcf_count(particleDistances[:,0], numOfBins, boxSize)
        P[j] = pressure(particleDistances)
        D[j] = diffusion(numOfParticles,numOfDimensions,initialX,PC3T[:,:,0,i])
        
    ########################
    ### Post-proecessing ###
    ########################

    if energyPlot == True:
        plotEnergy(numOfParticles,numOfTimesteps,timestep,saveFigures,U,T)
    #Pair correlation function
    pcf, pcfErr = pcf_plot(pcfCount, numOfParticles, numOfTimesteps, numOfBins, boxSize, saveFigures = False)
    # Pressure
    # P gives pressure at each timestep without time avg of prev timesteps
    P = (1 - 119.8/(3*numOfParticles*bathTemperature)*P)
    # PP gives mean pressure with time avg
    avgP, Perr = averageAndError(P,True, True)
    print('pressure is: {} $\ pm$ {}'.format(avgP, Perr))
    plt.figure('Pressure')
    xdata = np.arange(0,numOfTimesteps)
    plt.plot(xdata, P, color='b')
    plt.errorbar(xdata, np.full(xdata.shape,avgP), np.full(xdata.shape,Perr), color='r')
    plt.xlabel('t')
    plt.ylabel('P')
    plt.show()
    
    # Specific Heat (per particle)
    # The two formulas are derived from the same relation published by Verlet:
    # <dK^2>/<K> = T(1-3/(2C)) where C is the specific heat per particle
    
    # The lecture notes provide the formula:
    # <dK^2>/<K>^2 = 2/(3N)*(1-2/(2C))
    # which makes use of equipartition: <K> = 3/2 * N * T / 119.8
    
    # Since the mean kinetic energy will not always exactly match with our
    # user set temperature, there will be some discrepancy between the forms.
    
    # We can also use a 3rd form used by Lebowitz & Verlet in the same paper
    # <dK^2> = 3*N*T^2/2 * (1-3/(2C))
    # Here, the user set temperature is used instead of the average kinetic energy
    
    meanT, meanTErr = averageAndError(T,True,True)
    _, meanTsqErr = averageAndError(np.power(T,2)/1000,True,True) # the division by 1000 ensures the values don't get too large and curve_fit still works for finding the correlation length
    meanTsq = sum(np.power(T,2)) / len(T)
    varTsq = meanTsq - meanT**2
    #Cv1 = 1/((1 - varTsq/meanT**2 * 3*numOfParticles/2) * 2 / 3)
    #print("Specific heat (formula 1): ", Cv1)
    #Cv2 = 1/((1 - varTsq/meanT / (bathTemperature/119.8)) * 2 / 3)
    #print("Specific heat (formula 2): ", Cv2)
    Cv3 = 1/((1-varTsq*2/3/numOfParticles/(bathTemperature/119.8)**2) * 2 / 3)
    
    #print("meanTErr: ", meanTErr,", meanTsqErr: ", meanTsqErr)
    NTsq = numOfParticles*(bathTemperature/119.8)**2
    dCdKsq = (4 * 9 * NTsq) / (6 * NTsq - 4 * varTsq)**2
    dCdK = (-8 * meanT * 9 * NTsq) / (6 * NTsq - 4 * varTsq)**2
    dCv3 = np.sqrt((dCdKsq)**2 * meanTsqErr**2 + (dCdK)**2 * meanTErr**2)
    
    print('Specific heat: {} $\pm$ {}'.format(Cv3,dCv3))
    
    # use averageAndError(T,True,False) and averageAndError(np.power(T,2),True,False)
    # to compute errors and then add in quadrature
    
    # Diffusion
    plt.figure('Diffusion')
    time = np.arange(0, numOfTimesteps*timestep, timestep)
    plt.plot(time,D,'r')
    plt.xlabel('$3 \cdot 10^{-13}s$')
    plt.show()
    

    return MDS_dict,D
