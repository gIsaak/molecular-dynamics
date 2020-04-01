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

def getKineticEnergy(v):
    '''
    Accepts a velocity submatrix v containing all velocity components for all
    particles.

    The funciton returns the kinetic energy as calculated from the submatrix v.
    '''
    KE = np.sum(v**2)/2.0
    return KE

def getDistance(x,boxSize):
    '''
    Accepts a position submatrix x containing all position components for all n
    particles. Box size is also passed to enforce periodic boundary conditions

    Returns: nxn matrix Dist containing the distance d between each
            particle pair
            [d(1-1), d(2-1), ..., d(n-1)]
            [  ... ,  .... , ...,  .... ]
            [d(1-n),  .... , ..., d(n-n)]
            nx3xn matrix distComp containing the distance component between each
            particle pair
            (dim0 = particle1, dim1 = coordinate,  dim3 = particle2)
    '''
    expd = np.expand_dims(x,2)
    tiled = np.tile(expd, x.shape[0])
    trans = np.transpose(x)
    distComp = trans - tiled #exploitng numpy broadcasting
    distComp = np.mod((distComp + boxSize/2.0), boxSize) - boxSize/2.0 #periodic boundary conditions
    dist = np.sqrt(np.sum(distComp**2, axis=1))
    return dist, distComp


def getForce(dist,distComp):
    '''
    Function to get force on each particle for a given configuration

    Accepts dist(nxn) and distComp(nx3xn) arrays computed by getDistance
    The function returns a matrix F(nx3) containing the forces acting on the
    particles.
    '''
    N = dist.shape[0] #number of particles
    F = np.zeros((N,3))
    # get 1/r matrix
    invDist = np.reciprocal(dist, out=np.zeros_like(dist), where=dist!=0)
    # Computer force
    grad = 24*(invDist**8 - 2*invDist**14)
    grad_expand = np.expand_dims(grad,1)
    grad_tile = np.tile(grad_expand, (1,3,1))
    F = np.sum(np.multiply(distComp, grad_tile), axis=2) #broadcasting
    return F

def getPotentialEnergy(dist):
    '''
    Function to calculate total potential energy in the system.

    Accepts dist(nxn) array computed by getDistance
    The function returns a float PE of potential energy
    '''
    invDist = np.reciprocal(dist, out=np.zeros_like(dist), where=dist!=0)
    # Compute Potential Energy
    PE = 2*np.sum(invDist**12 - invDist**6)#1/2 to account for double counting of distances
    return PE

# Verlet algorithm
def iterateCoordinates_Verlet(x,v,f,timestep,boxSize):
    '''
    Accepts 3 submatrices for particle positions (x), particle velocities (v), and
    particle forces (f) as well as numbers of particles and dimensions and the size of the timestep
    and the box size. The function returns a new submatrix newX containing the Verlet-time-evolved
    positions at the next timestep.
    '''
    # To avoid redundant computations getForce is only called once per timestep in
    # iterateVelocities_Verlet
    newX = x + v*timestep + 0.5*f*timestep**2
    newX = np.mod(newX, boxSize) #periodic boundary condition
    return newX

# Verlet algorithm
def iterateVelocities_Verlet(v,f,nextf,timestep):
    '''
    Accepts 3 submatrices for particle velocities (v), particle forces at the same
    timestep as v (f), and particle forces at the timestep after v (nextf) as well
    as numbers of particles and dimensions and the size of the timestep
    The function returns a new submatrix newV containing the Verlet-time-evolved
    velocities at the next timestep.
    '''
    newV = v + 0.5*timestep*(f + nextf)
    return newV

def init_position(a,numOfParticles,density=0):
    '''
    Accepts lattice constant (a), numOfParticles
    Returns a submatrix newX of particles positions on fcc lattice
    as well as new number of particles (newNumOfParticles) and box size (L)

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

    newX = np.zeros((newNumOfParticles,3), dtype=float)
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

def gaussVel(temp, numOfParticles):
    '''
    Accepts temperature (temp), numOfParticles
    Returns submatrix newV of velocity components per particle initialized
    according to Gaussian distribution.
    The mean value of each velocity component array is set to zero.
    '''
    newV = np.zeros((numOfParticles,3), dtype=float)
    # Distribution parameters
    mu, sigma = 0, np.sqrt(temp/119.8) # mean is 0 and standard deviation in Kelvin
    for i in range(3):
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

def autocorr(data, plot=True, guess=[0,0]):
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

    guess: array
        supply [tau,a] guess for curve_fit routine
        this is important especially if dealing with large values, such
        as calculating the uncertainty of <T^2>

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
        a, b, c, d, e = 0, 0, 0, 0, 0
        #for n in range(N_sim - t):
        #    a = a + data[n]*data[n+t]
        #    b = b + data[n]
        #    c = c + data[n+t]
        #chi[t] = a/(N_sim - t) - b*c/(N_sim - t)**2
        for n in range(N_sim - t):
            a = a + data[n]*data[n+t]
            b = b + data[n]
            c = c + data[n+t]
            d = d + data[n]*data[n]
            e = e + data[n+t]*data[n+t]
        numChi = (N_sim - t)*a - b*c
        denChi = np.sqrt((N_sim - t)*d - b**2)*np.sqrt((N_sim - t)*e - c**2)
        chi[t] = numChi/denChi
    # Fit
    xdata = np.arange(chi_length)
    def func(x, tau, a):
        return np.exp(-x/tau)*a
    param = [0,0]
    if guess == [0,0]:
        param, _ = curve_fit(func, xdata, chi)
    else:
        param, _ = curve_fit(func, xdata, chi, guess)
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

def averageAndError(obs,autoCorr=False, plotAutoCorr=False, guess = [0,0]):
    '''
    Function to computer and observable average and error
    - autoCorr = False (default) assumes statistically independent data
    - autoCorr = True            rescales with correlation time tau

    Accepts: obs: array of observable values at each time-step
             numOfTimesteps: number of observable values
    Returns: avg: avergae of observable
             sigma: observable error
    '''
    obs2 = obs**2
    avg = np.mean(obs)
    avg2 = np.mean(obs2)
    N = obs.size
    sigma = np.sqrt((avg2 - avg**2)/(N - 1))
    if autoCorr == True:
        tau = 1
        if guess == [0,0]:
            _, tau = autocorr(obs, plotAutoCorr)
        else:
            _, tau = autocorr(obs,plotAutoCorr,guess)
        sigma = sigma*np.sqrt(2*tau*(N - 1)/N)
    return avg, sigma


def pcf_count(dist, numOfBins, boxSize, numOfParticles):
    '''
    Function to create histcount to compute PCF (Pair Correlation  Function)
    Accepts: dist: matrix from getDistance
             numOfBins: number of bins in hsitogram
             boxSize: simultion box size
    Returns: histCount array
             bin_edges array
    '''
    range = (0, boxSize)
    histCount, bin_edges = np.histogram(dist, bins=numOfBins, range=range)
    # take care of 0 entries in dist
    histCount[0] = histCount[0] - numOfParticles
    return histCount/2.0

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
    # plot only untill L/2
    half = np.size(r)//2
    plt.errorbar(r[:half], pcf[:half], pcfErr[:half], label='Pair Correlation Function')
    plt.xlabel('r')
    plt.ylabel('PCF')
    plt.grid()
    fig.suptitle('Pair correlation function', size = 14)
    # Save figure
    if saveFigures == True:
        plt.savefig('{}.png'.format(name), dpi=300)
    plt.show()
    return pcf, pcfErr


def pressure(dist):
    '''
    Function to compute pressure at time t
    Accepts: dist: matrix of particle distances
             numOfParticles: number of particles in the simulation
             bathTemperature: temperature
    Returns: P: pressure
    '''
    invDist = np.reciprocal(dist, out=np.zeros_like(dist), where=dist!=0)
    P = 6* np.sum(-2*invDist**12  + invDist**6) #1/2 factor due to double counting
    return P

def evolve(X, V, F, dist, distComp, timestep, boxSize):
    '''
    Function to let the system evolve one timestep
    Accepts: X, V, F, dist, distComp:
                    matrices of system parameters at current timestep
             timestep, boxSize:
                    simulation parameters
    Returns: Xnext, Vnext, Fnext, distNext, distCompNext:
                    matrices of system parameters at next timestep
    '''
    Xnext = iterateCoordinates_Verlet(X, V, F, timestep, boxSize)
    distNext, distCompNext = getDistance(Xnext, boxSize)
    Fnext = getForce(distNext, distCompNext)
    Vnext = iterateVelocities_Verlet(V, F, Fnext, timestep)
    return Xnext, Vnext, Fnext, distNext, distCompNext

def equilibrate(X,V,F,dist,distComp,timestep,boxSize,equilibrationTimer,bathTemperature):
    '''
    Function to set system at equilibrium

    Check if at equilibrium every equilibrationTimer steps
    If temperature agrees to desired temperature within <1> Kelvin, equilibrium is reached
    '''
    N = X.shape[0] #number of paricles
    equilibrium = False
    i = 0
    T_equilibration = 0
    while equilibrium == False:
        # Evolution
        X, V, F, dist, distComp = evolve(X, V, F, dist, distComp, timestep, boxSize)
        # Equilibrium control
        T = getKineticEnergy(V)
        # Equilibium check
        if i%equilibrationTimer == equilibrationTimer-1:
            if abs(float(T) / (N-1) / (3/2) * 119.8 - bathTemperature) > 1:
                scalingConstant = getVelocityScalingFactor(T, bathTemperature, N)
                print('Scaling constant:', scalingConstant)
                V = np.multiply(V, scalingConstant) #rescale
                T = getKineticEnergy(V) #update kinetic energy
                print('Rescaling from temperature: ',float(T) / (N-1) / (3/2) * 119.8)
            else:
                equilibrium = True
                print('Rescaled temperature: ',float(T) / (N-1) / (3/2) * 119.8)
                print('From initial temperature: {}'.format(bathTemperature))
        i += 1
    print('{} steps to equilibrium'.format(i))
    return X, V, F, dist, distComp

# TODO wait for ludi to fix this
#def diffusion(numOfParticles,initialX,PC3T):
#    '''
#    Function to gather position data of each particle to calculate diffusion
#    Accepts: numOfParticles: number of particles in the simulation
#             initialX: initial positions of particles
#             PC3T[:,:,0,i]: all particles x,y,z positions at timestep i
#    Returns: diff: array with mean square distance moved by all particles
#
#    '''
#    c_diff = np.zeros(shape=(numOfParticles,2))
#    diff = np.zeros(shape=(numOfParticles))
#    #for g in range(numOfParticles):
#    for d in range(3):
#        c_diff[:,0] = c_diff[:,0] + (PC3T[:,d])**2
#        c_diff[:,1] = c_diff[:,1] + (initialX[:,d])**2
#    diff = (np.sqrt(c_diff[:,0])-np.sqrt(c_diff[:,1]))**2
#    return np.mean(diff)

################# Begin main program ########################

def main(MDS_dict):

    #######################
    ### Load Dictionary ###
    #######################
    numOfTimesteps     = MDS_dict['number_of_timesteps']
    timestep           = MDS_dict['timestep']
    equilibrationTimer = MDS_dict['equilibration_timer']
    initParticles      = MDS_dict['init_particles']
    latticeConstant    = MDS_dict['lattice_constant']
    boxSize            = MDS_dict['box_size']
    numOfParticles     = MDS_dict['number_of_particles']
    bathTemperature    = MDS_dict['bath_temperature']
    plotting           = MDS_dict['plotting']
    plotCounter        = MDS_dict['plot_counter']
    energyPlot         = MDS_dict['energy_plot']
    density            = MDS_dict['density']
    saveFigures        = MDS_dict['save_figures']

    ########################
    ### Simulation Setup ###
    ########################
    # Initialize observables
    U = np.zeros((numOfTimesteps+1,1), dtype=float)
    T = np.zeros((numOfTimesteps+1,1), dtype=float)
    numOfBins = int(round(np.sqrt(numOfTimesteps))) #histogram rule
    pcfCount = np.zeros(shape=(numOfTimesteps,numOfBins))
    P = np.zeros(shape=(numOfTimesteps,))
    D = np.zeros(shape=(numOfTimesteps))
    D_avg = np.zeros(shape=(numOfTimesteps))
    # Miscellanea
    scatterPoints = []
    colours = ['b','g','r','c','m','y','k','w']

    # Initialization
    if initParticles == 'fcc':
        # place particles on fcc lattice and initialize velocities with M-B distribution
        X, numOfParticles, boxSize = init_position(latticeConstant, numOfParticles, density)
        V = gaussVel(bathTemperature, numOfParticles)
        dist, distComp  = getDistance(X, boxSize)
        F = getForce(dist, distComp)
    elif initParticles == 'random':
        raise Exception('Not coded yet')
    elif initParticles == 'debug':
        raise Exception('Not coded yet')
    # Plot
    if plotting == True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for p in range(numOfParticles):
            colour = colours[p%7]
            scatterPoints.append(ax.scatter(X[p,0],X[p,1],X[p,2],color=colour))
        ax.set_xlim((0,boxSize))
        ax.set_ylim((0,boxSize))
        ax.set_zlim((0,boxSize))
        plt.ion()
        plt.show()
        plt.pause(0.01)


    ##################
    ### Simulation ###
    ##################
    # Equilibration
    X, V, F, dist, distComp = \
            equilibrate(X, V, F, dist, distComp, timestep, boxSize, equilibrationTimer, bathTemperature)

    for i in range(numOfTimesteps):
        # Compute obsrvables
        U[i] = getPotentialEnergy(dist)
        T[i] = getKineticEnergy(V)
        pcfCount[i,:] = pcf_count(dist, numOfBins, boxSize, numOfParticles)
        P[i] = pressure(dist)
        #D[j] = diffusion(numOfParticles,initialX,PC3T[:,:,0,i])

        # Plot
        if plotting == True:
            if i%plotCounter == 0:
                for p in range(len(scatterPoints)):
                    scatterPoints.pop(0).remove()
                    colour = colours[p%7]
                    scatterPoints.append(ax.scatter(X[p,0],X[p,1],X[p,2],color=colour))
                plt.pause(0.000005)

        #Evolution
        X, V, F, dist, distComp = \
                evolve(X, V, F, dist, distComp, timestep, boxSize)

    ########################
    ### Post-proecessing ###
    ########################
    # Energy
    if energyPlot == True:
        plotEnergy(numOfParticles,numOfTimesteps,timestep,saveFigures,U,T)

    # Pair correlation function
    pcf, pcfErr = pcf_plot(pcfCount, numOfParticles, numOfTimesteps, numOfBins, boxSize, saveFigures)

    # Pressure
    # P gives pressure at each timestep without time avg of prev timesteps
    P = (1 - 119.8/(3*numOfParticles*bathTemperature)*P)
    # avgP gives mean pressure with time avg
    avgP, Perr = averageAndError(P,True, True)
    print('pressure is: {} $\ pm$ {}'.format(avgP, Perr))
    # Plot pressure
    plt.figure('Pressure')
    xdata = np.arange(0,numOfTimesteps)
    plt.plot(xdata, P, color='b')
    plt.errorbar(xdata, np.full(xdata.shape,avgP), np.full(xdata.shape,Perr), color='r')
    plt.xlabel('t')
    plt.ylabel('P')
    plt.show()

    # Heat Capacity (per particle)
    
    # We calculate the total heat capacity Cv according to the formula given in the lecture notes:
    # <dK^2>/<K>^2 = 2/(3N) * (1-3N/(2Cv))
    
    # This can be rearranged to f = <K^2>/<K>^2 = 1 + 2/(3N) - 1/Cv
    # we can use the known formula to compute the uncertainty of f and then propagate
    # this uncertainty to Cv
    
    meanT, meanTErr = averageAndError(T,True,True,[10,1000])
    meanTsq, meanTsqErr = averageAndError(np.power(T,2),True,True,[10,1.5e9])
    
    f = meanTsq/meanT**2
    covKKsq = np.mean(np.multiply((T - meanT),(np.power(T,2)-meanTsq)))
    df = abs(f)*np.sqrt((meanTsqErr/meanTsq)**2 + (2*meanTErr/meanT)**2 - 2*covKKsq/meanTsq/meanT**2)
    
    Cv = (1+2/(3*numOfParticles)-f)**(-1)
    dCv = df/(1+2/(3*numOfParticles)-f)**2
    
    # convert to heat capacity per particle
    C = Cv/numOfParticles
    dC = dCv/numOfParticles
    print('Heat Capacity: {} +/- {}'.format(C,dC)))

    # Diffusion
    #plt.figure('Diffusion')
    #time = np.arange(0, numOfTimesteps*timestep, timestep)
    #plt.plot(time,D,'r')
    #plt.xlabel('$3 \cdot 10^{-13}s$')
    #plt.show()

    return MDS_dict
