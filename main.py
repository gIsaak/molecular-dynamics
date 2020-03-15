import simulation_func as fu



# initialisation options:
#   1) 'fcc'        --YET TO IMPLEMENT
#   2) 'random'     --hanldes any number of particles
#   3) 'debug_2'    --requires the numOfParticles to be 2

# plot_counter=int  --integer value to display only every int's scatter plot to reduce simulation plots

# savefig == True   --saves figure as MDS_#particles_method.png in current dir
'''
Parameters

algorithm:              {euler, verlet}


latticeConstant:        integer, fcc lattice constant in units of sigma
numOfParticles :        integer, 

'''


MDS_parameters = {

# Simulation
'algorithm':            'verlet',
'number_of_timesteps':  1000,
'timestep' :            0.001,

# Initialial Configuration
'init_particles':       'fcc',
'box_size':             6, # init_particles 'fcc' will overwrite
'lattice_constant':     3, # used only by init_particles 'fcc'
'number_of_particles':  14, # init_particles 'fcc' will overwrite
'number_of_dimensions': 3, # hardcoded to 3 in simulation_func.py force calculation
'bath_temperature' :    300, # temperature of microcanonical ensemble

# Plotting options
'plotting':             True, # Plot particle motion in scatter plot
'plot_counter' :        10,
'energy_plot' :         True,
'save_figures':         False,

# Observables
'pair_correlation':     True
# TODO create other observable flags
}

_ = fu.main(MDS_parameters)