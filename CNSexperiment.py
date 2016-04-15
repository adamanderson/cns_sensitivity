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
#   --all distances are in cm

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
import scipy.integrate as spint
import pdb
import ReactorTools

# fundamental constants
keVPerGeV       = 1e6      # [keV / GeV]
keVPerMeV       = 1e3
gPerkg          = 1e3
hbarc 			= 0.197*keVPerGeV		# [keV fm]
fmPercm			= 1.0e13	# [fm / cm]
sin2thetaW		= 0.2387
cLight			= 3.0e8 	# [m / s]
nAvogadro		= 6.022e23
joulePereV		= 1.602e-19	# [J / eV]
eVPerFission	= 200.0e6 	# [eV]
Mn              = 0.931 * keVPerGeV
Gfermi			= (1.16637e-5 / (keVPerGeV**2.))*(hbarc/fmPercm)		# [cm / keV]


class CNSexperiment:
    def __init__(self, N, Z, dRdEnu, dRdT_b, detMass, time, nuFlux):
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
        self.nNeutrons 		    = N
        self.nProtons 		    = Z
        self.nNucleons          = N + Z
        self.isotopeMassc2	    = self.nNucleons * Mn 		# [keV]
        self.dRdEnu_source      = dRdEnu
        self.dRdT_background    = dRdT_b
        self.Qweak				= self.nNeutrons - (1.0 - 4.0*sin2thetaW) * self.nProtons
        self.detectorMass       = detMass
        self.nTargetAtoms       = self.detectorMass * gPerkg / self.nNucleons * nAvogadro # approximately
        self.livetime           = time
        self.nuFlux             = nuFlux


    def dsigmadT_atEnu_CNS(self, Enu, T):
        '''

        Parameters
        ----------


        Returns
        -------
        dRdE : float
            Differential rate at given energy
        '''
        if type(Enu) == float:
    		Enu = np.asarray([Enu])
    	else:
    		Enu = np.asarray(Enu)
        if type(Enu) == float:
    		Enu = np.asarray([Enu])
    	else:
    		Enu = np.asarray(Enu)
        dsigmadT = Gfermi**2. / (4*np.pi) * self.Qweak**2. * self.isotopeMassc2 * \
                    (1. - (self.isotopeMassc2 * T) / (2.*Enu**2.)) * self.F_Helm(T, self.nNucleons)
        dsigmadT[dsigmadT<0] = 0.
        return dsigmadT


    def dsigmadT_CNS(self, T):
        '''
        docs
        '''
        def dsigmadTdEnu_CNS(Enu, T):
            dsigmadTdEnu = 1.e42 * self.dsigmadT_atEnu_CNS(Enu, T) * self.dRdEnu_source(Enu)
            # print 'dsigmadT = %f' % self.dsigmadT_atEnu_CNS(Enu, T)
            # print 'nu spectrum = %f' % self.dRdEnu_source(Enu)
            return dsigmadTdEnu

        dsigmadT = np.array([spint.quad(dsigmadTdEnu_CNS, 0, 1.e6, args=(Tval), epsabs=1.e-6)[0]/1.e42 for Tval in T])
        return dsigmadT


    def dRdT_CNS(self, T):
        dRdT = self.dsigmadT_CNS(T) * self.nuFlux * self.livetime * self.nTargetAtoms
        return dRdT

    def F_Helm(self, T, A):
        # define the momentum transfer in MeV / c
        q = np.sqrt(2 * A * Mn * T)

        # Standard Helm form factor
        R = 0.89*A**(1.0/3.0) + 0.3   # nuclear radius [fm]
        a = 0.52                      # [fm]
        s = 0.9                       # smearing parameter [fm]
        c = 1.23 * A**(1.0/3.0) - 0.6 # [fm]
        rn = np.sqrt(c**2.0 + (7.0/3.0)*(np.pi*a)**2 - 5*s**2)   # [fm]
        F = (3 / ((q * rn / hbarc)**3) * (np.sin(q * rn / hbarc) - (q * rn / hbarc) * np.cos(q * rn / hbarc)))**2.0 * \
            np.exp(-1.0 * (q*s / hbarc)**2)

        return F


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
