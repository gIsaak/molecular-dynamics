import simulation_func as fu



# initialisation options:
#   1) 'fcc'        --YET TO IMPLEMENT
#   2) 'random'     --hanldes any number of particles
#   3) 'debug_2'    --requires the numOfParticles to be 2

# plot_counter=int  --integer value to display only every int's scatter plot to reduce simulation plots

# savefig == True   --saves figure as MDS_#particles_method.png in current dir



MDS_parameters = {
'euler' :           False,
'verlet' :          True,
'boxSize_L' :       5,
'numOfParticles' :  2,
'numOfDimensions' : 3,
'num_t' :           20000,
'timestep' :        0.001,
'plotting':         True,
'plot_counter' :    50,
'energyPlot' :      True,
'save_fig' :        False,
'init_particles' :  'debug_2'
}


_ = fu.main(MDS_parameters)

