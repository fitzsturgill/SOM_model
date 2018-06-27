
from pylab import *
import utils

N_Py = 25
N_SOM = 1
N_PV = 1
stim_rf = array([0])



class circuit():
    # eventualy I think I want to separate the weights from some of the things that are more plot oriented- i.e. receptive field
    # the circuit is the most basic thing
    # how the circuit is then run and plotted are 2 additional things
    # i.e. you might create a circuit instance and modify the weights
    # You then might want to plot it across different ranges of stimulus strength
    # for optimal code reuse would you use the same instance but modify the stimulus strength ranges?

    def __init__(self):
## create circuit instance
        # Py = Excitatory, additive
        self.N_Py = N_Py
        self.Py_rf = linspace(-180, 180, N_Py)

        # SOM = Inhibitory, additive
        self.N_SOM = N_SOM

        # PV = Inhibitory, multiplicative
        self.N_PV = N_PV

        ## Baseline rates -- light off
        self.bPy = 0 #1
        self.bSOM = 0 #1
        self.bPV = 0 #1

        # Synaptic weights
        self.wInPy = 0.5 * utils.cos_weight(stim_rf, Py_rf, 1)
        self.wInSOM = 0. * ones((1, N_SOM))
        self.wInPV = 0. * ones((1, N_PV))

        self.wPyPy = 0. * ones((N_Py, N_Py))
        self.wPySOM = 5. / N_Py * ones((N_Py, N_SOM))
        self.wPyPV = 5. / N_Py * ones((N_Py, N_PV)) #5. / N_Py * ones((N_Py, N_PV))

        self.wSOMPy = -2. * ones((N_SOM, N_Py)) # -2.
        self.wSOMSOM = -0. * ones((N_SOM, N_SOM))
        self.wSOMPV = -0.1 * ones((N_SOM, N_PV)) #-0.05 * ones((N_SOM, N_PV))

        # default is for PV cells to be divisive
        self.additive = 0

    def run(self, rateIn):

        """
            calc_equilibrium: calculates the equilibrium firing rate for a feedback network
            of recti-linear rate neurons.
            """


        ratePy = zeros((self.N_Py, 1))
        rateSOM = zeros((self.N_SOM, 1))
        ratePV = zeros((self.N_PV, 1))

        numsteps = 1000
        stepsize = .01

        #ratePy_hist = zeros(numsteps)
        #rateSOM_hist = zeros(numsteps)

        for i in range(numsteps):

            if self.additive:
                ratePy_eq = (rateIn * self.wInPy.T + dot(self.rateSOM.T, self.wSOMPy) + dot(self.ratePy.T, self.wPyPy) + self.bPy + dot(self.ratePV.T, self.wPVPy))
            else:
                ratePy_eq = (rateIn * self.wInPy.T + dot(self.rateSOM.T, self.wSOMPy) + dot(self.ratePy.T, self.wPyPy) + self.bPy)  / (dot(self.ratePV.T, self.wPVPy) + 1)

            rateSOM_eq = (rateIn * self.wInSOM.T + dot(self.ratePy.T, self.wPySOM) + dot(self.rateSOM.T, self.wSOMSOM) + self.bSOM) / (dot(self.ratePV.T, self.wPVSOM) + 1)

            ratePV_eq = (rateIn * self.wInPV.T + dot(self.ratePy.T, self.wPyPV) + dot(self.rateSOM.T, self.wSOMPV) + bPV) / (dot(self.ratePV.T, self.wPVPV) + 1)

            ratePy += stepsize * (ratePy_eq.T - ratePy)
            rateSOM += stepsize * (rateSOM_eq.T - rateSOM)
            ratePV += stepsize * (ratePV_eq.T - ratePV)

            ratePy[ratePy < 0] = 0
            rateSOM[rateSOM < 0] = 0
            ratePV[ratePV < 0] = 0

        return ratePy, rateSOM, ratePV        


class rateModel():

    def __init__(self, circuit, stim_strength):
        self.PyRates = zeros((circuit.N_Py, len(stim_strength), len(circuit.stim_rf)))
        self.SOMRates = zeros((circuit.N_SOM, len(stim_strength), len(circuit.stim_rf)))
        self.PVRates = zeros((circuit.N_PV, len(stim_strength), len(circuit.stim_rf)))


        ## Run the simulation

        for i in arange(len(stim_strength)):
            rateIn = stim_strength[i]        
            ratePy_eq, rateSOM_eq, ratePV_eq = circuit.run(rateIn)            
            for j in arange(len(circuit.stim_rf)):
                self.PyRates[:, i, j] = ratePy_eq[:, 0]
                self.SOMRates[:, i, j] = rateSOM_eq[:, 0]
                self.PVRates[:, i, j] = ratePV_eq[:, 0]

