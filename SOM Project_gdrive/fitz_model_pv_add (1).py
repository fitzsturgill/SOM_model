# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 14:44:28 2014

@author: Paxon
"""


from pylab import *
import utils


def calc_equilibrium(rateIn, bPy, bSOM, bPV, wInPy, wInSOM, wInPV, wPyPy, wPySOM, wPyPV, wSOMPy, wSOMSOM, wSOMPV, wPVPy, wPVSOM, wPVPV, additive):
    """
    calc_equilibrium: calculates the equilibrium firing rate for a feedback network
    of recti-linear rate neurons.
    """
    
    N_Py = wInPy.shape[0]
    N_SOM = wInSOM.shape[1]
    N_PV = wInPV.shape[1]
    
    ratePy = zeros((N_Py, 1))
    rateSOM = zeros((N_SOM, 1))
    ratePV = zeros((N_PV, 1))
    
    numsteps = 1000
    stepsize = .01
    
    #ratePy_hist = zeros(numsteps)
    #rateSOM_hist = zeros(numsteps)
    
    for i in range(numsteps):
        
        #ratePy_eq = (dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy))  / (dot(ratePV.T, wPVPy) + 1)
        ratePy_eq = (rateIn * wInPy.T + dot(rateSOM.T, wSOMPy) + dot(ratePy.T, wPyPy) + bPy + dot(ratePV.T, wPVPy))
        if i==1 and additive:
            print("hello")
        
        rateSOM_eq = (rateIn * wInSOM.T + dot(ratePy.T, wPySOM) + dot(rateSOM.T, wSOMSOM) + bSOM + dot(ratePV.T, wPVSOM))
        
        ratePV_eq = (rateIn * wInPV.T + dot(ratePy.T, wPyPV) + dot(rateSOM.T, wSOMPV) + bPV + dot(ratePV.T, wPVPV))

        ratePy += stepsize * (ratePy_eq.T - ratePy)
        rateSOM += stepsize * (rateSOM_eq.T - rateSOM)
        ratePV += stepsize * (ratePV_eq.T - ratePV)
        
        ratePy[ratePy < 0] = 0
        rateSOM[rateSOM < 0] = 0
        ratePV[ratePV < 0] = 0
                            
    return ratePy, rateSOM, ratePV
    

###############################################################################


## Set up the circuitry

# The stimulus type -- receptive field
stim_rf = array([0])

# The stimulus strength
stim_strength = arange(20, 101, 5)

# Py = Excitatory, additive
N_Py = 25
Py_rf = linspace(-180, 180, N_Py)

# SOM = Inhibitory, additive
N_SOM = 1

# PV = Inhibitory, multiplicative
N_PV = 1

# Baseline rates -- light off
bPy = 0
bSOM = 0
bPV = 0

# Synaptic weights -- light off
wInPy = 0.5 * utils.cos_weight(stim_rf, Py_rf, 1)
wInSOM = 0. * ones((1, N_SOM))
wInPV = 0. * ones((1, N_PV))

wPyPy = 0. * ones((N_Py, N_Py))
wPySOM = 5. / N_Py * ones((N_Py, N_SOM))
wPyPV = 5. / N_Py * ones((N_Py, N_PV))

wSOMPy = -2 * ones((N_SOM, N_Py))
wSOMSOM = -0. * ones((N_SOM, N_SOM))
wSOMPV = -0.05 * ones((N_SOM, N_PV))

wPVPy = -1.0 * ones((N_PV, N_Py))
wPVSOM = 0. * ones((N_PV, N_SOM))
wPVPV = 0. * ones((N_PV, N_PV))

# Baseline rates -- light on
bPy0 = 0.
bSOM0 = 0.
bPV0 = 0.

# Synaptic weights -- light on
wInPy0 = 1. * utils.cos_weight(stim_rf, Py_rf, 1)
wInSOM0 = 0. * ones((1, N_SOM))
wInPV0 = 0. * ones((1, N_PV))

wPyPy0 = 0. * ones((N_Py, N_Py))
wPySOM0 = 1. / N_Py * ones((N_Py, N_SOM))
wPyPV0 = 5. / N_Py * ones((N_Py, N_PV))

wSOMPy0 = -1.5 * ones((N_SOM, N_Py))
wSOMSOM0 = -0. * ones((N_SOM, N_SOM))
wSOMPV0 = -0.5 * ones((N_SOM, N_PV))

wPVPy0 = -1.0* ones((N_PV, N_Py))
wPVSOM0 = 0. * ones((N_PV, N_SOM))
wPVPV0 = 0. * ones((N_PV, N_PV))


ratePy_eqs = zeros((N_Py, len(stim_strength), len(stim_rf)))
rateSOM_eqs = zeros((N_SOM, len(stim_strength), len(stim_rf)))
ratePV_eqs = zeros((N_PV, len(stim_strength), len(stim_rf)))

ratePy0_eqs = zeros((N_Py, len(stim_strength), len(stim_rf)))
rateSOM0_eqs = zeros((N_SOM, len(stim_strength), len(stim_rf)))
ratePV0_eqs = zeros((N_PV, len(stim_strength), len(stim_rf)))


## Run the simulation

for i in arange(len(stim_strength)):
    for j in arange(len(stim_rf)):
        
        rateIn = stim_strength[i]
        
        ratePy_eq, rateSOM_eq, ratePV_eq = calc_equilibrium(rateIn, bPy, bSOM, bPV,
                            wInPy, wInSOM, wInPV, wPyPy, wPySOM, wPyPV, 
                            wSOMPy, wSOMSOM, wSOMPV, wPVPy, wPVSOM, wPVPV, 1)
                            
        
        ratePy0_eq, rateSOM0_eq, ratePV0_eq = calc_equilibrium(rateIn, bPy0, bSOM0, bPV0,
                            wInPy0, wInSOM0, wInPV0, wPyPy0, wPySOM0, wPyPV0, 
                            wSOMPy0, wSOMSOM0, wSOMPV0, wPVPy0, wPVSOM0, wPVPV0, 1)
                            
                        
        ratePy_eqs[:, i, j] = ratePy_eq[:, 0]
        rateSOM_eqs[:, i, j] = rateSOM_eq[:, 0]
        ratePV_eqs[:, i, j] = ratePV_eq[:, 0]
                     
        ratePy0_eqs[:, i, j] = ratePy0_eq[:, 0]
        rateSOM0_eqs[:, i, j] = rateSOM0_eq[:, 0]
        ratePV0_eqs[:, i, j] = ratePV0_eq[:, 0]
        
        
## Save the data
        
today = datetime.date.today()
        
filename = 'fitz_model-pv_add-' + today.strftime('%y%m%d')       

savez(filename, stim_strength=stim_strength, stim_rf=stim_rf, 
      ratePy_eqs=ratePy_eqs, rateSOM_eqs=rateSOM_eqs, ratePV_eqs=ratePV_eqs,
      ratePy0_eqs=ratePy0_eqs, rateSOM0_eqs=rateSOM0_eqs, ratePV0_eqs=ratePV0_eqs,
      bPy=bPy, bSOM=bSOM, bPV=bPV, wInPy=wInPy, wInSOM=wInSOM, wInPV=wInPV,
      wPyPy=wPyPy, wPySOM=wPySOM, wPyPV=wPyPV, wSOMPy=wSOMPy, wSOMSOM=wSOMSOM, wSOMPV=wSOMPV,
      wPVPy=wPVPy, wPVSOM=wPVSOM, wPVPV=wPVPV,
      bPy0=bPy0, bSOM0=bSOM0, bPV0=bPV0, wInPy0=wInPy0, wInSOM0=wInSOM0, wInPV0=wInPV0,
      wPyPy0=wPyPy0, wPySOM0=wPySOM0, wPyPV0=wPyPV0, wSOMPy0=wSOMPy0, wSOMSOM0=wSOMSOM0, wSOMPV0=wSOMPV0,
      wPVPy0=wPVPy0, wPVSOM0=wPVSOM0, wPVPV0=wPVPV0)
         
        
## Load the data
        
#np.load('fitz_model_noPV-141013.npz')
     
## Plot the data

c_idxs = array([12]) # which Py cells to plot
#c_idxs = array([12])
#a_idxs = array([0, 1, 2]) # which touch angles to plot
a_idxs = arange(len(stim_rf))
i_idxs = array([len(stim_strength)-3]) # which intensities to plot


matplotlib.rcParams.update({'font.size':20})

# Plot pyramidal cells
figure(1)
clf()

for a in a_idxs:
    for c in c_idxs:

        plot(stim_strength, squeeze(ratePy_eqs[c, :, a]), c='k', lw=5, label='LED Off')
        plot(stim_strength, squeeze(ratePy0_eqs[c, :, a]), c=[1, 0.2, 0.2], lw=5, label='LED On')


ymin, ymax = ylim()

xlabel('Intensity/Concentration (arb)')
ylabel('Firing Rate (spikes/s)')
title('Pyramidal')

legend(loc='upper left')

f1_filename = filename + '-stim_py'
savefig(f1_filename + '.png', bbox_inches='tight')
savefig(f1_filename + '.eps', bbox_inches='tight')


# Plot inhibitory cells                     
f2 = figure(2)
clf()

subplot(211)
plot(stim_strength, squeeze(rateSOM_eqs[0, :, 0]), c='k', lw=4, label='LED Off')
plot(stim_strength, squeeze(rateSOM0_eqs[0, :, 0]), c=[1, 0.2, 0.2], lw=4, label='LED On')

legend(loc='upper left')
title('SOM')
#ylabel('Firing Rate (spikes/s)')

subplot(212)
plot(stim_strength, squeeze(ratePV_eqs[0, :, 0]), c='k', lw=4, label='PV')
plot(stim_strength, squeeze(ratePV0_eqs[0, :, 0]), c=[1, 0.2, 0.2], lw=4)

#legend(loc='upper left')
title('PV')
xlabel('Intensity/Concentration (arb)')
#ylabel('Firing Rate (spikes/s)')
f2.text(0.06, 0.5, 'Firing Rate (spikes/s)', ha='center', va='center', rotation='vertical')


f2_filename = filename + '-stim_SOMPV'
savefig(f2_filename + '.png', bbox_inches='tight')
savefig(f2_filename + '.eps', bbox_inches='tight')

# Plot pyramidal cell receptive fields
figure(3)
clf()

#rcPynarams.update({'font.size': 20})

for a in a_idxs:
    for i in i_idxs:
        
        plot(Py_rf, squeeze(ratePy_eqs[:, i, a]), c='k', lw=5, label='LED Off')
        plot(Py_rf, squeeze(ratePy0_eqs[:, i, a]), c=[1, 0.2, 0.2], lw=5, label='LED On')
        

ylim((ymin, ymax))
xlim((-180, 180))

xlabel('Receptive Field')
ylabel('Firing Rate (spikes/s)')
title('Pyramidal')

legend(loc='upper left')

f3_filename = filename + '-rf_py'
savefig(f3_filename + '.png', bbox_inches='tight')
savefig(f3_filename + '.eps', bbox_inches='tight')

# Plot Py cell firing against each other
figure(4)
clf()

plot(squeeze(ratePy_eqs[c_idxs[0], :, a_idxs[0]]),
     squeeze(ratePy0_eqs[c_idxs[0], :, a_idxs[0]]), 'ok', ms=5)
plot([ymin, ymax], [ymin, ymax], c=[0.7, 0.7, 0.7])

axis('scaled')
xlabel('LED Off (spikes/s)')
ylabel('LED On (spikes/s)')
title('Py Firing Rate')


f4_filename = filename + '-py_py'
savefig(f4_filename + '.png', bbox_inches='tight')
savefig(f4_filename + '.eps', bbox_inches='tight')




