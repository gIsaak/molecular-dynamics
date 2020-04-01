import simulation_func_vector as fu
import diffusion_processor as diff
import h5py

##### User-populated parameters/arrays for batch simulations #####
T_085= [0.591,0.658,0.719,0.760,0.782,0.786,0.880,1.128,1.214,2.202,2.889,3.000,3.500,4.174]
T_045 = [1.55, 1.71, 2.93, 4.62]
T_010 = [5.5,6.5,7.5]
rho = 0.45

# Initialisation of arrays for diffusion calculation
# must be individually adjusted to the loop calling the main function
ite = 20 # number of iterations to average diffusion curve
d_saver_10 = np.zeros(shape=(10000)) # shape = number_of_timesteps
diffErr = np.zeros(shape=(ite))

C = []
dC = []
P = []
dP = []


for t in T_045: 
    MDS_parameters = {
    
    # Simulation
    'number_of_timesteps':  10000,   # Number of evolutions performed after equilibrium is reached.
    'timestep' :            0.001,   # Size of timestep. Note unity in natural units equals 2.15e-12 physical seconds
    'equilibration_timer':  100,     # Number of simultion steps before velocity rescaling during equilibration.
    
    # Initialial Configuration
    'init_particles':       'fcc',   # Only to be modified for debugging purposes
    'box_size':             6,       # init_particles 'fcc' will overwrite. Only used for debugging purposes
    'lattice_constant':     1.49,    # Used only by init_particles 'fcc' when density set to 0
    'number_of_particles':  365,     # init_particles 'fcc' will overwrite by rounding to: 14,63,172,365,etc. to fill lattice
    'bath_temperature' :    119.8*t, # Temperature of microcanonical ensemble, use t as dimensionless temp
    'density':			    rho,     # Used only by init_particles 'fcc'. Choose 0 for arbitrary density set by user specification of a
    
    # Plotting options
    'plotting':             False,   # Plot particle motion in scatter plot
    'plot_counter' :        10,      # Number of evolutions before scatter plot updated
    'energy_plot' :         True,    # Plot of kinetic, potential, and total energies after equilibrium
    'save_figures':         False,   # Flag to automatically save generated figures
    
    }
    avgP,Perr,avgC,Cerr,D,Derr = fu.main(MDS_parameters)
    C.append(avgC)
    dC.append(Cerr)
    P.append(avgP)
    dP.append(Perr)
    d_saver_10 += D
    diffErr[t] = Derr

# optional saving of data in binary file
#hf = h5py.File('simulation_add_parameter_info.h5', 'w')
#hf.create_dataset('spec_heat', data=C)
#hf.create_dataset('spec_heat_err', data=dc)
#hf.create_dataset('pressure', data=P)
#hf.create_dataset('pressure_err', data=dp)
#hf.create_dataset('diffusion', data=d_saver_10)
#hf.create_dataset('diffusion_err', data=diffErr)
#hf.close()

# funtion to plot diffusion for 3 differnt diffusion data sets
# it is adviseable to work in diffusion_processor seperately to generate the plots
#values = diff.plot_diffusion(G,GErr,L,LErr,S,SErr)
