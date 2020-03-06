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
'boxSize_L' :       5, #optional, init_position will overwirte
'latticeConst' :    3,
'numOfParticles' :  12,
'numOfDimensions' : 3,
'temp' :            300,
'num_t' :           500,
'timestep' :        0.001,
'plotting':         True,
'plot_counter' :    1,
'energyPlot' :      True,
'save_fig' :        False,
'init_particles' :  'fcc'
}


_ = fu.main(MDS_parameters)
