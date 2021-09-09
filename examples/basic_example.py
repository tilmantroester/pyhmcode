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
mass_def = ccl.halos.MassDef("vir", 'matter')
hmf = ccl.halos.MassFuncSheth99(ccl_cosmology, mass_def=mass_def,
                                mass_def_strict=False, use_delta_c_fit=True)
hbf = ccl.halos.HaloBiasSheth99(ccl_cosmology, mass_def=mass_def,
                                mass_def_strict=False, use_delta_c_fit=True)
hmc = ccl.halos.HMCalculator(ccl_cosmology, hmf, hbf, mass_def)

# Compute the CCL halo model power spectrum, using the halo profile from HMx
ccl_halomodel_pofk = ccl.halos.halomod_Pk2D(
                            cosmo=ccl_cosmology,
                            hmc=hmc, 
                            prof=profile_generator.matter_profile,
                            normprof1=False,
                            a_arr=a, lk_arr=np.log(k*h))