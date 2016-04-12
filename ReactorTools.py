# ReactorTools.py
#
# Some tools for calculating the neutrino rate from a nuclear reactor.
#
# Adam Anderson
# 14 April 2016
# adama@fnal.gov
#
# Note: Convention on units:
#   --all masses are in kg
#   --all energies are in keV


def dRdEnu_U235(Enu):
	'''
	Reactor anti neutrino spectrum from U235 (see arXiv:1101.2663v3)

	Parameters
	----------
	Enu : array
		Neutrino energy in MeV

	Returns
	-------
	spectrum : array
		Spectrum [nu / MeV / fission]
	'''
	spectrum = np.exp(3.217 - 3.111*Enu + 1.395*(Enu**2.0) - \
					  (3.690e-1)*(Enu**3.0) + (4.445e-2)*(Enu**4.0) - (2.053e-3)*(Enu**5.0))
	spectrum[Enu<1.0] = np.exp(3.217 - 3.111*1.0 + 1.395*(1.0**2.0) - \
					  (3.690e-1)*(1.0**3.0) + (4.445e-2)*(1.0**4.0) - (2.053e-3)*(1.0**5.0))
	return spectrum


def dRdEnu_U238(Enu):
	'''
	Reactor anti neutrino spectrum from U238 (see arXiv:1101.2663v3)

	Parameters
	----------
	Enu : array
		Neutrino energy in MeV

	Returns
	-------
	spectrum : array
		Spectrum [nu / MeV / fission]
	'''
	spectrum = np.exp((4.833e-1) + (1.927e-1)*Enu - (1.283e-1)*Enu**2.0 - \
						(6.762e-3)*Enu**3.0 + (2.233e-3)*Enu**4.0 - (1.536e-4)*Enu**5.0)
	spectrum[Enu<1.0] = np.exp((4.833e-1) + (1.927e-1)*1.0 - (1.283e-1)*1.0**2.0 - \
						(6.762e-3)*1.0**3.0 + (2.233e-3)*1.0**4.0 - (1.536e-4)*1.0**5.0)
	return spectrum


def dRdEnu_Pu239(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3)

	Parameters
	----------
	Enu : array
		Neutrino energy in MeV

	Returns
	-------
	spectrum : array
		Spectrum [nu / MeV / fission]
	'''
	spectrum = np.exp(6.413 - 7.432*Enu + 3.535*Enu**2.0 - \
						(8.82e-1)*Enu**3.0 + (1.025e-1)*Enu**4.0 - (4.550e-3)*Enu**5.0)
	spectrum[Enu<1.0] = np.exp(6.413 - 7.432*1.0 + 3.535*1.0**2.0 - \
						(8.82e-1)*1.0**3.0 + (1.025e-1)*1.0**4.0 - (4.550e-3)*1.0**5.0)
	return spectrum
