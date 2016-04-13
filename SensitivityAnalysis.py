import CNSexperiment
import ReactorTools
import matplotlib.pyplot as plt
import numpy as np

# setup
def bkg(T): return 0
myExpt = CNSexperiment.CNSexperiment(N=40, Z=32, dRdEnu=ReactorTools.dRdEnu_U235, dRdT_b=bkg)

T = np.linspace(0, 10, 100)
dsigmadT_1MeV = myExpt.dsigmadT_atEnu_CNS(1000., T)

# plots
Enu = np.linspace(0, 1e4, 100)
plt.figure()
plt.plot(Enu, ReactorTools.dRdEnu_U235(Enu), label='U235')
plt.plot(Enu, ReactorTools.dRdEnu_U238(Enu), label='U238')
plt.plot(Enu, ReactorTools.dRdEnu_Pu239(Enu), label='Pu239')
plt.legend()
plt.xlabel('energy [keV]')
plt.ylabel('nu / keV / fission')

plt.figure()
plt.plot(T, dsigmadT_1MeV)
plt.xlabel('recoil energy [keV]')
plt.ylabel('')

plt.show()
