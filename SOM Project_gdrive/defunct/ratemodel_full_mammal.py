# -*- coding: utf-8 -*-
"""
Created on Mon May 28 21:39:50 2012

This is a script to test out whether the feedback loop can lead to a stable
state with the multiplicative gain. This is simplified to a rate model to get
a better idea of just the feedback dynamics.

@author: Paxon
"""

from brian import *
import utils

def calc_equilibrium(rateIn, wInPy, wInSOM, wPySOM, wSOMPy, plothist=False):
    """
    Calculates the equilibrium rates of the local bend In -> interneuron system.
    """
    
    N_Py = wInPy.shape[1]
    N_SOM = wInSOM.shape[1]
    
    ratePy = zeros((N_Py, 1))
    rateSOM = zeros((N_SOM, 1))
    
    ratePysum = 0.0
    InPysum = 0.0
    
    #eps = 1e-20
        
    numsteps = 1000
    stepsize = .01
    
    #ratePy_hist = zeros(numsteps)
    #rateSOM_hist = zeros(numsteps)
    
    for i in range(numsteps):
        
        ratePysum = 0.0
        
        # Additive synapse
        ratePy_eq = dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy)        
        # Multiplicative synapse
        #rateI_eq = dot(rateIn.T, wInI) / (dot(rateSOM.T, wSOMI) + 1)
        #rateI_eq[k] = float(wInI[k] * rateIn) / (wSOMI[k] * rateSOM + 1)
        ##


        rateSOM_eq = dot(rateIn.T, wInSOM) + dot(ratePy.T, wPySOM)

        ratePy += stepsize * (ratePy_eq.T - ratePy)
        rateSOM += stepsize * (rateSOM_eq.T - rateSOM)
        
        ratePy[ratePy < 0] = 0
        rateSOM[rateSOM < 0] = 0
        
                    
        ## SOMABA Saturates
        #rateSOM[rateSOM > 1] = 1            
        
    return ratePy, rateSOM
       
def calc_equilibrium_x(rateIn, wInPy, wInSOM, wPySOM, wSOMPy, wInPV, wPyPV, wPVPy, wSOMPV=0, wPVSOM=0, SOMbaseline=0, plothist=False):
    """
    Calculates the equilibrium rates of the local bend In -> interneuron system.
    """
    
    N_Py = wInPy.shape[1]
    N_SOM = wInSOM.shape[1]
    N_PV = wInPV.shape[1]
    
    ratePy = zeros((N_Py, 1))
    rateSOM = zeros((N_SOM, 1))
    ratePV = zeros((N_PV, 1))
    
    ratePysum = 0.0
    InPysum = 0.0
    
    #eps = 1e-20
        
    numsteps = 1000
    stepsize = .01
    
    #ratePy_hist = zeros(numsteps)
    #rateSOM_hist = zeros(numsteps)
    
    for i in range(numsteps):
        
        ratePysum = 0.0
        
        #ratePy_eq = (dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy))  / (dot(ratePV.T, wPVPy) + 1)
        ratePy_eq = (dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy))  / (dot(ratePV.T, wPVPy) + 1)
        rateSOM_eq = dot(rateIn.T, wInSOM) + dot(ratePy.T, wPySOM) + dot(ratePV.T, wPVSOM) +  SOMbaseline
        
        ratePV_eq = dot(rateIn.T, wInPV) + dot(ratePy.T, wPyPV) + dot(rateSOM.T, wSOMPV) + 3

        ratePy += stepsize * (ratePy_eq.T - ratePy)
        rateSOM += stepsize * (rateSOM_eq.T - rateSOM)
        ratePV += stepsize * (ratePV_eq.T - ratePV)
        
        ratePy[ratePy < 0] = 0
        rateSOM[rateSOM < 0] = 0
        ratePV[ratePV < 0] = 0
                            
        ## SOMABA Saturates
        #rateSOM[rateSOM > 1] = 1            
        
    return ratePy, rateSOM, ratePV
       
       

def calc_equilibrium_0(rateIn, wInPy, wInSOM, wPySOM, wSOMPy, wInPV, wPyPV, wPVPy, wSOMPV=0, wPVSOM=0, SOMbaseline=0, plothist=False):
    """
    Calculates the equilibrium rates of the local bend In -> interneuron system.
    """
    
    N_Py = wInPy.shape[1]
    N_SOM = wInSOM.shape[1]
    N_PV = wInPV.shape[1]
    
    ratePy = zeros((N_Py, 1))
    rateSOM = zeros((N_SOM, 1))
    ratePV = zeros((N_PV, 1))
    
    ratePysum = 0.0
    InPysum = 0.0
    
    #eps = 1e-20
        
    numsteps = 1000
    stepsize = .01
    
    #ratePy_hist = zeros(numsteps)
    #rateSOM_hist = zeros(numsteps)
    
    for i in range(numsteps):
        
        ratePysum = 0.0
        
        #ratePy_eq = (dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy))  / (dot(ratePV.T, wPVPy) + 1)

        ratePy_eq = (dot(rateIn.T, wInPy) + dot(rateSOM.T, wSOMPy))  / (dot(ratePV.T, wPVPy) + 1)

        rateSOM_eq =  0 * (dot(rateIn.T, wInSOM) + dot(ratePy.T, wPySOM) + dot(ratePV.T, wPVSOM)) + SOMbaseline
        
        ratePV_eq = dot(rateIn.T, wInPV) + dot(ratePy.T, wPyPV) + dot(rateSOM.T, wSOMPV) + 3

        ratePy += stepsize * (ratePy_eq.T - ratePy)
        rateSOM += stepsize * (rateSOM_eq.T - rateSOM)
        ratePV += stepsize * (ratePV_eq.T - ratePV)
        
        ratePy[ratePy < 0] = 0
        rateSOM[rateSOM < 0] = 0
        ratePV[ratePV < 0] = 0
                            
        ## SOMABA Saturates
        #rateSOM[rateSOM > 1] = 1            
        
    return ratePy, rateSOM, ratePV
              
def sigmoid(x):
  return 1 / (1 + exp(-0.1*x + 5))              
              
    ################################################################################

touch_angle = array([0]) # This can only handle 1 angle for now.
#touch_angle = arange(0, 91, 15)
touch_strength = arange(0, 51, 2)

N_In = 4
In_rf = [-135, -45, 45, 135]

N_Py = 25
#Py_rf = arange(-180.0, 180.0, 360.0/N_Py)
Py_rf = linspace(-180, 180, N_Py)

#N_SOM = 4
N_SOM = 1
N_PV = 1


InPy_weight = 10
wInPy = InPy_weight * utils.cos_weight(Py_rf, In_rf, 1)
#wInI = InI_weight * utils.lin_weight(I_rf, In_rf, 1)

InSOM_weight = 0
wInSOM = InSOM_weight * ones((N_In, N_SOM))
#wInSOM = InSOM_weight * eye(N_In, N_SOM)
wPySOM = 0.5 * ones((N_Py, N_SOM))
# Uniform SOM I weight distribution
wSOMPy = -2 * ones((N_SOM, N_Py))
# Uniform for mult gain
#wSOMPy = 0 * ones((N_SOM, N_Py))
# Opposite cos distribution
#wSOMPy = -0.7 * utils.cos_weight(Py_rf, Pyn_rf, 1)

wInPV = 0 * ones((N_In, N_PV))

wPyPV = 5.0 / (N_Py) * ones((N_Py, N_PV))

wPVPy = 0.1 * ones((N_PV, N_Py))

wPVSOM = 0

wSOMPV = -0.05 * ones((N_SOM, N_PV))

## Shift the excitatory weight matrix to see what happens with a different cell
#wInI = append(wInI[3:], wInI[0:3])
##
rateIn_eqs = zeros((N_In, len(touch_strength), len(touch_angle)))
ratePy_eqs = zeros((N_Py, len(touch_strength), len(touch_angle)))
rateSOM_eqs = zeros((N_SOM, len(touch_strength), len(touch_angle)))
ratePV_eqs = zeros((N_SOM, len(touch_strength), len(touch_angle)))

ratePy0_eqs = zeros((N_Py, len(touch_strength), len(touch_angle)))
rateSOM0_eqs = zeros((N_SOM, len(touch_strength), len(touch_angle)))
ratePV0_eqs = zeros((N_SOM, len(touch_strength), len(touch_angle)))

for i in arange(len(touch_strength)):    
    for j in arange(len(touch_angle)):
        #rateIn = touch_strength[i] * utils.cos_weight(touch_angle, In_rf, 1)   
        rateIn = touch_strength[i] * utils.cos_weight(array([touch_angle[j]]), In_rf, 1)    

        #rateIn = touch_strength[i] * utils.lin_weight(touch_angle, In_rf, 2)
    
        #rateI_eq, rateSOM_eq = calc_equilibrium(rateIn, wInI, wInSOM, wISOM, wSOMI)
        #rateI0_eq, rateSOM0_eq = calc_equilibrium(rateIn, wInI, wInSOM, wISOM, (wSOMI*0))  
        ratePy_eq, rateSOM_eq, ratePV_eq = calc_equilibrium_x(rateIn, wInPy, wInSOM, wPySOM, wSOMPy, wInPV, wPyPV, wPVPy, wSOMPV, wPVSOM, 20)
        ratePy0_eq, rateSOM0_eq, ratePV0_eq = calc_equilibrium_0(rateIn, wInPy, wInSOM, wPySOM, (wSOMPy), wInPV, wPyPV, (wPVPy), (wSOMPV), wPVSOM, 0)
        
        rateIn_eqs[:, i, j] = rateIn[:, 0]
        ratePy_eqs[:, i, j] = ratePy_eq[:, 0]
        rateSOM_eqs[:, i, j] = rateSOM_eq[:, 0]
        ratePV_eqs[:, i, j] = ratePV_eq[:, 0]
        
        ratePy0_eqs[:, i, j] = ratePy0_eq[:, 0]
        rateSOM0_eqs[:, i, j] = rateSOM0_eq[:, 0]
        ratePV0_eqs[:,i , j] = ratePV0_eq[:, 0]


# Inlot parameters
c_idxs = array([12]) # which cells to plot
#c_idxs = array([12])
#a_idxs = array([0, 1, 2]) # which touch angles to plot
a_idxs = arange(len(touch_angle))
i_idxs = array([8]) # which intensities to plot

figure(3)
clf()

#a = gca()
#a.set_xticks([])
#a.set_yticks([])

#rcParams.update({'font.size': 20})


for a in a_idxs:
    for c in c_idxs:
        col = utils.angle2color(Py_rf[c], 90)
        col = col[0,:]
        #col = 'k'
#        
#        plot(touch_strength, squeeze(rateI_eqs[c, :, a]), c=col, lw=5, label='Normal')
#        plot(touch_strength, squeeze(rateI0_eqs[c, :, a]), '--', c=col, lw=5, label='SOMABA Blocked')

        plot(touch_strength, squeeze(ratePy_eqs[c, :, a]), c=col, lw=5, label='Cell %0i$\degree$' % Py_rf[c] )
        plot(touch_strength, squeeze(ratePy0_eqs[c, :, a]), '--', c=col, lw=5)
        
        for i in i_idxs:
            plot(touch_strength[i], ratePy_eqs[c, i, a], 'o', c=col, ms=15)
            plot(touch_strength[i], ratePy0_eqs[c, i, a], '^', c=col, ms=15)


figure(5)
clf();
plot(touch_strength, squeeze(rateSOM_eqs[0, :, 0]), c=[0, 0, 0.5], lw=4, label='SOM')
plot(touch_strength, squeeze(rateSOM0_eqs[0, :, 0]), '--', c=[0, 0, 0.5], lw=4)

plot(touch_strength, squeeze(ratePV_eqs[0, :, 0]), c=[0.5, 0, 0], lw=4, label='PV')
plot(touch_strength, squeeze(ratePV0_eqs[0, :, 0]), '--', c=[0.5, 0, 0], lw=4)

legend(loc='upper left')
ymin, ymax = ylim()

figure(6)
clf()

a = gca()
a.set_xticks([])
a.set_yticks([])

#rcPynarams.update({'font.size': 20})

for a in a_idxs:
    for i in i_idxs:
        col = utils.angle2color(touch_angle[a], 90)
        col = col[0, :]
        
        plot(Py_rf, squeeze(ratePy_eqs[:, i, a]), c=col, lw=5,)
        plot(Py_rf, squeeze(ratePy0_eqs[:, i, a]), '--', c=col, lw=5)
        
        for c in c_idxs:
            col = utils.angle2color(Py_rf[c], 90)
            col = col[0,:]
            plot(Py_rf[c], ratePy_eqs[c, i, a], 'o', c=col, ms=15)
            plot(Py_rf[c], ratePy0_eqs[c, i, a], '^', c=col, ms=15)

ylim((ymin, ymax))
