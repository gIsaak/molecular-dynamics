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
'boxSize_L' :       6, #optional, init_position will overwirte
'latticeConst' :    1.7,
'numOfParticles' :  14,
'numOfDimensions' : 3,
'temp' :            90,
'num_t' :           5000,
'timestep' :        0.001,
'plotting':         True,
'plot_counter' :    10,
'energyPlot' :      True,
'save_fig' :        False,
'init_particles' :  'fcc'
}


_ = fu.main(MDS_parameters)


# When I test with 14 particles in a box of size 50, energy is conserved just fine

# When debug 3 is run, there are sharp spikes in the kinetic and total energies.
# Aside from these individual spikes, total energy remains quite constant.

# BIG PROBLEM FOUND: In getForce function, potential energy was overwritten for each
# pairwise interaction, not summed over the total number of pairwise interaction