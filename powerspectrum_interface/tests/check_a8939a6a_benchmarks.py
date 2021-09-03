# Check against the version of the library and interface used for the paper
# (commit a8939a6a)

import os
import sys
import numpy as np

import camb
import pyhmx


if __name__ == "__main__":
    create_benchmarks = len(sys.argv) > 1 and sys.argv[1] == "create"

    hmx = pyhmx.HMx()

    # Cosmological parameters for CAMB
    h = 0.7
    omc = 0.25
    omb = 0.048
    mnu = 0.06
    w = -1.0
    wa = 0.0
    ns = 0.97
    As = 2.1e-9

    z_max = 2.0
    n_z = 8

    k_max = 20.0

    # Get linear power spectrum
    # Set up CAMB
    p = camb.CAMBparams(WantTransfer=True,
                        WantCls=False,
                        Want_CMB_lensing=False,
                        DoLensing=False,
                        NonLinearModel=None,
                        )
    p.set_cosmology(H0=h*100, omch2=omc*h**2, ombh2=omb*h**2, mnu=mnu)
    p.set_dark_energy(w=w, wa=wa)
    p.set_initial_power(camb.InitialPowerLaw(As=As, ns=ns))

    z_lin = np.linspace(0, z_max, n_z, endpoint=True)
    p.set_matter_power(redshifts=z_lin, kmax=k_max, nonlinear=False)

    r = camb.get_results(p)

    # Get sigma_8, linear power spectrum, and Omega_m
    sigma8 = r.get_sigma8()[-1]
    k_lin, z_lin, pofk_lin_camb = r.get_matter_power_spectrum(
                                                minkh=1e-3,
                                                maxkh=20.0,
                                                npoints=128)

    omv = r.omega_de + r.get_Omega("photon") + r.get_Omega("neutrino")
    omm = p.omegam

    cosmology = {"Omega_m"  : omm,
                 "Omega_b"  : omb,
                 "Omega_v"  : omv,
                 "h"        : h,
                 "n_s"      : ns,
                 "sigma_8"  : sigma8,
                 "m_nu"     : mnu,
                 "w"        : w,
                 "wa"       : wa}

    halo_model = {"eta0" : 0.603,
                  "As"   : 3.13,
                  "Theat": 10**7.8}

    halo_model_mode = pyhmx.constants.HMx2020_matter_pressure_with_temperature_scaling

    # Run HMCode
    Pk_HMx_dmonly = hmx.run_HMCode(cosmology=cosmology,
                                   halo_model=halo_model,
                                   k=k_lin,
                                   z=z_lin,
                                   pk_lin=pofk_lin_camb,
                                   verbose=True)
    
    # Run HMCode
    Pk_HMx = hmx.run_HMx(cosmology=cosmology,
                         halo_model=halo_model,
                         fields=[pyhmx.constants.field_matter,
                                 pyhmx.constants.field_electron_pressure],
                         mode=halo_model_mode,
                         k=k_lin,
                         z=z_lin,
                         pk_lin=pofk_lin_camb,
                         verbose=True)
    
    output_dir = "a8939a6a_benchmarks"

    if not create_benchmarks:
        z_lin_bm = np.loadtxt(os.path.join(output_dir, "pofk_lin_z.txt"))
        k_lin_bm = np.loadtxt(os.path.join(output_dir, "pofk_lin_k.txt"))
        pofk_lin_camb_bm = np.loadtxt(os.path.join(output_dir, "pofk_lin.txt"))
        Pk_HMx_dmonly_bm = np.loadtxt(os.path.join(output_dir, "pofk_hmcode2016.txt"))
        Pk_HMx_mm_bm = np.loadtxt(os.path.join(output_dir, "pofk_hmx_matter_matter.txt"))
        Pk_HMx_mp_bm = np.loadtxt(os.path.join(output_dir, "pofk_hmx_matter_pressure.txt"))

        assert np.allclose(Pk_HMx_dmonly, Pk_HMx_dmonly_bm, rtol=2e-3)
        assert np.allclose(Pk_HMx_mm_bm, Pk_HMx[0, 0], rtol=2e-3)
        assert np.allclose(Pk_HMx_mp_bm, Pk_HMx[0, 1], rtol=6e-3)
    else:
        os.makedirs(output_dir, exist_ok=True)
        np.savetxt(os.path.join(output_dir, "pofk_lin_z.txt"), z_lin)
        np.savetxt(os.path.join(output_dir, "pofk_lin_k.txt"), k_lin)
        np.savetxt(os.path.join(output_dir, "pofk_lin.txt"), pofk_lin_camb)
        np.savetxt(os.path.join(output_dir, "pofk_hmcode2016.txt"), Pk_HMx_dmonly)
        np.savetxt(os.path.join(output_dir, "pofk_hmx_matter_matter.txt"), Pk_HMx[0,0])
        np.savetxt(os.path.join(output_dir, "pofk_hmx_matter_pressure.txt"), Pk_HMx[0,1])

