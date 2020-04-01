import numpy as np
import h5py
import matplotlib.pyplot as plt
  
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def SI_to_dimless(numOfParticles, SI_density):
    '''
    Function to get dimensionless density from g/cm^3 density
    accepts: numOfParticles: number of particles
             SI_density: Density in SI units
    
    returns: Density in dimensionless units
    
    '''
    m_Ar = 39.95*1.6474*10**-24
    x_cm3 = numOfParticles*m_Ar/SI_density
    x_A3 = x_cm3*10**24
    L = (x_A3)**(1/3)/3.405
    
    return numOfParticles/(L)**3
    

def plot_diffusion(G,GErr,L,LErr,S,SErr,numOfTimesteps=20000,ts=0.001,ite=20):
    '''
    Function to plot the diffusion for three different diffusion ismulations 
    (denoted as G-gas, L-liquid, S-solid)
    accepts: G: data array of gas simulations (not averaged yet)
             GErr: data array of errors of gas simulation (not averaged yet)
             L: data array of liquid simulations (not averaged yet)
             LErr: data array of errors of liquid simulation (not averaged yet)
             S: data array of solid simulations (not averaged yet)
             SErr: data array of errors of solid simulation (not averaged yet)
             
             numOfTimesteps: number of timesteps used, default at 20000
             ts: time increment step size, default 0.001
             ite: number of simulations to average over
    
    returns: Diffusion coefficients in cm^2/s
    '''
    
    #linear fit
    time = np.arange(0, numOfTimesteps*ts, ts)
    G_diff,G_cov = np.polyfit(time[1500:4000],G[1500:4000]*(3.405)**2/ite,1,cov=True)
    L_diff,L_cov = np.polyfit(time[13000:],L[13000:]*(3.405)**2/ite,1,cov=True)
    S_diff,S_cov = np.polyfit(time[13000:],S[13000:]*(3.405)**2/ite,1,cov=True)
    
    #linear fit coefficients
    a = G_diff[1]
    G_diff = G_diff[0]
    
    b = L_diff[1]
    L_diff = L_diff[0]
    
    c = S_diff[1]
    S_diff = S_diff[0]
    
    #calculation of mean error plus linear fit error
    G_meanErr = np.sqrt(np.sum(GErr*(3.405**2)**2))/ite + G_cov[1,1]
    L_meanErr = np.sqrt(np.sum(LErr*(3.405**2)**2))/ite + L_cov[1,1]
    S_meanErr = np.sqrt(np.sum(SErr*(3.405**2)**2))/ite + S_cov[1,1]
    print('Fitting error Liquid: {}'.format(L_cov[1,1]))
    print('Fitting error Solid: {}'.format(S_cov[1,1]))
    
    # plotting
    
    ####### Diffusion - Gas with inset #########
    fig, ax = plt.subplots()
    ax.plot(time[:],G[:]*(3.405)**2/ite,label=r'T=300K $\rho$=0.01')
    plt.xlabel('time unit in $10^{-12}s$')
    plt.ylabel('$10^{-16}cm^2$')
    plt.legend()
    # first subplot
    plt.title('Diffusion - Gas')
    axins = inset_axes(ax,
                       width="60%", 
                       height="60%",
                       bbox_to_anchor=(0,-0.28,1,1), bbox_transform=ax.transAxes)
    axins.plot(time[:1250],G[:1250]*(3.405)**2/ite)
    axins.plot([0,time[1250]],[60,60],color='black',label='linear growth threshold')
    axins.legend()
    axins.grid()
    plt.xticks(visible=True)
    plt.yticks(visible=True)
    plt.legend()
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    plt.draw()
#    plt.savefig('G_diff_20avg_inset.png', dpi=500)
    plt.show()
    ####### end plot with inset #########
    
    ####### plot with linear fit #######
#    plt.figure('Diffusion')
#
#    plt.title('Diffusion - Liquid, Solid')
#    
#    #data and fit gas
##    plt.plot(time[:],G[:]*(3.405)**2/ite, label='mean of {} simulations N=365'.format(ite))
##    plt.plot(time[500:6000],G_diff*time[500:6000]+a, label='D = {:1.2f}$\pm${:1.2f}E-3'.format(G_diff/6/10,G_meanErr/6))
#    
#    # data and fit liquid
#    plt.plot(time[:],L[:]*(3.405)**2/ite, color='#ff7f0eff', label=r'T=94.4K $\rho$=0.81')
#    plt.plot(time[3000:],L_diff*time[3000:]+b, label='$D_{liquid}$'+' = {:1.2f}$\pm${:1.2f}$\cdot e-5$'.format(L_diff/6*(10),L_meanErr/6))
#    
#    # data and fit solid
#    plt.plot(time[:],S[:]*(3.405)**2/ite, color='#2ca02cff', label=r'T=10K $\rho$=0.81')
#    plt.plot(time[3000:],S_diff*time[3000:]+c, color='red', label='$D_{solid}$'+' = {:1.2f}$\pm${:1.2f}$\cdot e-5$'.format(S_diff/6*10,S_meanErr/6))
#    
#    plt.xlabel('time unit in $10^{-12}s$')
#    plt.ylabel('$10^{-16}cm^2$')
#    plt.legend()
#    plt.grid()
#    plt.savefig('LS_diff_20avg.png', dpi=500)
#    plt.show()
    ####### end plot with linear fit #########
    
    return [[G_diff/6, G_meanErr/6],[L_diff/6,L_meanErr/6],[S_diff/6,S_meanErr/6]]
    
## usage example 

#Load Data 
#data = h5py.File(r'w_dir\name.h5', 'r')
#G = data.get('diffusion')
#GErr = data.get('diffusion_err')
##...
#
#values = plot_diffusion(G,GErr,L,LErr,S,SErr)




