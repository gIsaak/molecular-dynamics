import simulation_func_vector as fu



# initialisation options:
#   1) 'fcc'        --YET TO IMPLEMENT
#   2) 'random'     --hanldes any number of particles
#   3) 'debug_2'    --requires the numOfParticles to be 2

# plot_counter=int  --integer value to display only every int's scatter plot to reduce simulation plots

# savefig == True   --saves figure as MDS_#particles_method.png in current dir

i = 0.1 # dimless temp as used in verlet paper

MDS_parameters = {

# Simulation
'number_of_timesteps':  5000,
'timestep' :            0.001,
'equilibration_timer':  20,      # number of simultion steps before velocity rescaling

# Initialial Configuration
'init_particles':       'fcc',
'box_size':             6,       # init_particles 'fcc' will overwrite
'lattice_constant':     1.49,    # used only by init_particles 'fcc'
'number_of_particles':  120,      # init_particles 'fcc' will overwrite
'bath_temperature' :    119.8*i, # temperature of microcanonical ensemble use i as dimless temp
'density':				0.06,    # used only by init_particles 'fcc' choose 0 for arbitrary density set by user specification of a

# Plotting options
'plotting':             False,    # Plot particle motion in scatter plot
'plot_counter' :        10,
'energy_plot' :         True,
'save_figures':         False,

}
_ = fu.main(MDS_parameters)
