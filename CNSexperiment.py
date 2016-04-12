# CNSexperiment.py
#
# Class that represents a CNS experiment with some tools for generating fake
# data and computing limits.
#
# Adam Anderson
# 12 April 2016
# adama@fnal.gov
#
# Note: Convention on units:
#   --all masses are in kg
#   --all energies are in keV

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import pdb

# fundamental constants
hbarc 			= 0.197		# [GeV fm]
fmPercm			= 1.0e13	# [fm / cm]
sin2thetaW		= 0.2387
cLight			= 3.0e8 	# [m / s]
nAvogadro		= 6.022e23
joulePereV		= 1.602e-19	# [J / eV]
eVPerFission	= 200.0e6 	# [eV]
massMuon		= 105.658   # [MeV]

class CNSexperiment:
    # particle parameters
    nNeutrons 		= 40.
    nProtons 		= 32.
    molarMass 		= 78.	# [g / mol]
    isotopeMassc2	= molarMass * 0.9314941 		# [GeV]
    Qw				= nNeutrons - (1.0 - 4.0*sin2thetaW) * nProtons
    Gfermi			= (1.16637e-5)*hbarc/fmPercm		# [cm / GeV]

    # detector parameters
    detectorMass 	= 1 * 1e3 	# [g]
    nNuclei 		= nAvogadro * (detectorMass/molarMass)

    def __init__(self, dRdEnu_signal, dRdT_background):
        '''
        Parameters
        ----------
        dRdEnu_signal : function
            Differential energy spectrum for signal neutrinos (e.g. reactor
            neutrino spectrum), **per neutrino energy**
        dRdT_background : function
            Differential energy spectrum for background, **per recoil energy**

        Returns
        -------
        None
        '''

    def dRdT_CNS(self, T):
        '''
        Differential spectrum as a function of energy for CNS.

        Parameters
        ----------
        T : float
            Energy in keV at which to evaluate the spectrum

        Returns
        -------
        dRdE : float
            Differential rate at given energy
        '''
        diffCS = np.zeros((len(T), len(Enu)))
    	for jEnu in range(len(Enu)):
    		if np.sum(T < Enu[jEnu]) > 0 and Enu[jEnu]>0:
    			diffCS[T < Enu[jEnu], jEnu] = (2.495e-25) * (kappa**2) * (1.0/T[T < Enu[jEnu]] - 1.0/Enu[jEnu])
    	return diffCS


    def dEdT_background(self, T):
        '''
        Differential spectrum as a function of energy for CNS.

        Parameters
        ----------
        T : float
            Energy in keV at which to evaluate the spectrum

        Returns
        -------
        dRdE : float
            Differential rate at given energy
        '''


    def neg2logL(self, mu):
        '''
        The -2logL as a function of a signal strength parameter mu. The
        parameter mu is defined such that mu=1 corresponds to the nominal cross
        section of the CNS signal.

        Parameters
        ----------
        mu : float
            Signal strength parameter

        Returns
        -------
        n2LogL : float
        '''


    def run_toy(self):
        '''
        Generates a pseudoexperiment by toy MC.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''


    def UL_test_stat(self, data):
        '''
        Test statistic for an upper limit.

        Parameters
        ----------
        data : array
            Numpy array with list of event energies

        Returns
        -------
        q0 : float
            Test statistic for upper limit
        '''
