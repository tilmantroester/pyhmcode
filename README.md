# pyhmcode

This library provides a convenient interface to Alexander Mead's library, which underlies the HMCode and HMx codes to predict non-linear power spectra. For details on the Fortran version, refer to https://github.com/alexander-mead/HMcode.
The `pyhmcode` interface uses the excellent [f90wrap](https://github.com/jameskermode/f90wrap) library to generate the python interface. This allows access to virtually all functionality in the library, although to speed up compilation, the interface is limited to commonly-used subroutines.
If you use `pyhmcode`, please cite [Tröster, Mead et al. 2021](https://arxiv.org/abs/2109.04458) and the relevant papers describing the model:
- HMCode2015: Mead et al. (2015; https://arxiv.org/abs/1505.07833)
- HMCode2016: Mead et al. (2016; https://arxiv.org/abs/1602.02154)
- HMCode2020: Mead et al. (2021; https://arxiv.org/abs/2009.01858)
- HMx: Mead, Tröster et al. (2020; https://arxiv.org/abs/2005.00009)

## Installation
`pyhmcode` is `pip`-installable, so
```
pip install pyhmcode
```
should do the trick. For more control over the installation and access to the alternative interface, as well as CosmoSIS support, the repository can be cloned with
```
git clone --recursive https://github.com/tilmantroester/pyhmcode
```

## Alternative interface and CosmoSIS support
There is an alternative interface to compute powerspectra in the `powerspectrum_interface` subdirectory. This interface uses `ctypes` and a thin Fortran wrapper to access the main functionality of HMCode and HMx, namely predicting the non-linear power spectra.
Installation proceeds by
```
cd powerspectrum_interface
pip install .
```
This installs the `pyhmx` python module, which in turn is interfaced with CosmoSIS using the `cosmosis_interface.py` module.

## Demos

The `notebooks` and `example` subdirectories include a number of examples on how to use the the interfaces. 
A basic example (from `examples/basic_example.py`) showing the use of `pyhmcode` to generate non-linear power specta, as well as its integration with CCL is shown below.
```python
import pyhmcode
import pyhmcode.halo_profile_utils

# We use CCL to generate the linear power spectrum
import pyccl as ccl

import numpy as np


ccl_cosmology = ccl.CosmologyVanillaLCDM()

k = np.logspace(-4, 1.5, 100)
a = np.linspace(1/(1+6), 1, 10)
z = 1/a - 1

pofk_lin = np.array([ccl.linear_matter_power(ccl_cosmology, k=k, a=a_)
                     for a_ in a])

# CCL uses units of Mpc, while pyhmcode uses Mpc/h. Hence we need to convert
# the units here.
h = ccl_cosmology["h"]
k = k/h
pofk_lin = pofk_lin * h**3

# Create the pyhmcode cosmology object. Beside the cosmology parameters, it
# also holds the linear power spectrum.
hmcode_cosmology = pyhmcode.halo_profile_utils.ccl2hmcode_cosmo(
                        ccl_cosmo=ccl_cosmology,
                        pofk_lin_k_h=k,
                        pofk_lin_z=z[::-1],
                        pofk_lin=pofk_lin[::-1],
                        log10_T_heat=7.8)

# Create the halo model object, which holds information on the specific halo
# model to use. E.g., the HMCode or HMx version.
hmcode_model = pyhmcode.Halomodel(
                    pyhmcode.HMx2020_matter_pressure_w_temp_scaling)

# Now we can compute the non-linear power spectrum, given the cosmology,
# halo model, and a list of fields.
hmcode_pofk = pyhmcode.calculate_nonlinear_power_spectrum(
                                    cosmology=hmcode_cosmology,
                                    halomodel=hmcode_model, 
                                    fields=[pyhmcode.field_matter,
                                            pyhmcode.field_electron_pressure])

# The output of calculate_nonlinear_power_spectrum has
# shape (n_field, n_field, n_z, n_k).
matter_matter_pofk = hmcode_pofk[0, 0]
matter_electron_pressure_pofk = hmcode_pofk[0, 1]

# We can also use the halo profiles from HMCode or HMx and use them in another
# halo model code.
profile_generator = pyhmcode.halo_profile_utils.HMxProfileGenerator(
                        hmcode_cosmology,
                        a_arr=a, k_arr=k,
                        fields=[pyhmcode.field_matter,
                                pyhmcode.field_cdm,
                                pyhmcode.field_electron_pressure],
                        add_diffuse=False)

# Here we use the halo profile in the CCL halo model framework. 
# First setup the CCL halo model specification
mass_def = ccl.halos.MassDef(Delta="vir", rho_type="matter")
hmf = ccl.halos.MassFuncSheth99(mass_def=mass_def,
                                mass_def_strict=False, use_delta_c_fit=True)
hbf = ccl.halos.HaloBiasSheth99(mass_def=mass_def,
                                mass_def_strict=False, use_delta_c_fit=True)
hmc = ccl.halos.HMCalculator(halo_bias=hbf, mass_function=hmf, mass_def=mass_def)

# Compute the CCL halo model power spectrum, using the halo profile from HMx
ccl_halomodel_pofk = ccl.halos.halomod_Pk2D(
                            cosmo=ccl_cosmology,
                            hmc=hmc, 
                            prof=profile_generator.matter_profile,
                            a_arr=a, lk_arr=np.log(k*h))
```
