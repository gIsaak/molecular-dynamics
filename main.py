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
'latticeConst' :    2,
'numOfParticles' :  60,
'numOfDimensions' : 3,
<<<<<<< HEAD
'num_t' :           500,
'timestep' :        0.001,
'plotting':         False,
'plot_counter' :    1,
=======
'num_t' :           20,
'timestep' :        0.001,
'plotting':         False,
'plot_counter' :    50,
>>>>>>> e06eb502672b3269b82048ab79006f7b94b79c7e
'energyPlot' :      True,
'save_fig' :        False,
'init_particles' :  'fcc'
}


_ = fu.main(MDS_parameters)
