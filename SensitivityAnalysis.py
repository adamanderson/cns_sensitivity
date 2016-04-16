import CNSexperiment
import ReactorTools
import matplotlib.pyplot as plt
import numpy as np

# setup
bkg = np.zeros(100)
T_bkg = np.linspace(0,100,100)
flux = ReactorTools.nuFlux(5.5, 400.)
myExpt = CNSexperiment.CNSexperiment(N=14, Z=14, dRdEnu=ReactorTools.dRdEnu_U235,
                                     T_bg=T_bkg, dRdT_bg=bkg, detMass=5, time=365*24*60*60., nuFlux=flux)

T = np.linspace(1e-6, 0.5, 100)
Thigh = np.linspace(1e-6, 100, 100)
dsigmadT_1MeV = myExpt.dsigmadT_atEnu_CNS(1000., T)
dsigmadT_U235 = myExpt.dsigmadT_CNS(T)
dRdT_U235 = myExpt.dRdT_CNS(T)
toydata = myExpt.run_toy()

mu = np.linspace(0.9, 1.1, 20)
neg2LogL = np.array([myExpt.neg2logL(m, toydata) for m in mu])

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

plt.figure()
plt.semilogy(T, dsigmadT_U235)
plt.xlabel('recoil energy [keV]')
plt.title('spectrum from U235 reactor neutrinos')
plt.ylabel('')

plt.figure()
plt.hist(toydata, bins=np.linspace(0, 0.5, 51))
plt.plot(T, dRdT_U235 / 100)
plt.xlabel('recoil energy [keV]')
plt.title('rate from U235 reactor neutrinos')
plt.ylabel('events / 10 eV')

plt.figure()
plt.plot(mu, neg2LogL)
plt.xlabel('signal / SM')
plt.ylabel('-2logL')
plt.show()
