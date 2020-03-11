import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def normal_autocorr(mu, sigma, tau, N):
    """Generates an autocorrelated sequence of Gaussian random numbers.

    Each of the random numbers in the sequence of length `N` is distributed
    according to a Gaussian with mean `mu` and standard deviation `sigma` (just
    as in `numpy.random.normal`, with `loc=mu` and `scale=sigma`). Subsequent
    random numbers are correlated such that the autocorrelation function
    is on average `exp(-n/tau)` where `n` is the distance between random
    numbers in the sequence.

    This function implements the algorithm described in
    https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf

    Parameters
    ----------

    mu: float
        mean of each Gaussian random number
    sigma: float
        standard deviation of each Gaussian random number
    tau: float
        autocorrelation time
    N: int
        number of desired random numbers

    Returns:
    --------
    sequence: numpy array
        array of autocorrelated random numbers
    """
    f = np.exp(-1./tau)

    sequence = np.zeros(shape=(N,))

    sequence[0] = np.random.normal(0, 1)
    for i in range(1, N):
        sequence[i] = f * sequence[i-1] + np.sqrt(1 - f**2) * np.random.normal(0, 1)

    return mu + sigma * sequence

def autocorr(data, plot=True):
    '''
    Computes autocorrelation function for an observable

    See formula in the lecutre notes.
    The cutoff for t is given by t = sqrt(N_sim)
    Parameters
    ----------

    data: numpy array
        value of observable at each timestep

    Returns:
    --------
    chi: numpy array
        autocorrelation function
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
    def func(x, tau, b):
        return np.exp(-x/tau) + b
    param, _ = curve_fit(func, xdata, chi)
    tau = param[0]
    print(tau)
    if plot == True:
        plt.plot(xdata, chi, 'b-', label='data')
        # Fit
        plt.plot(xdata, func(xdata, *param), 'r--', label='fit: tau=%5.3f, b=%5.3f' % tuple(param))
        plt.xlabel('t')
        plt.ylabel('Chi')
        plt.legend()
        plt.show()


    return chi, tau

N = 20000
mu = 0
sigma = 1
tau = 50
data = normal_autocorr(mu, sigma, tau, N)
chi = autocorr(data)

plt.plot(data)
plt.show()
