import numpy as np

import camb
import pyhmx

import contextlib
import time


@contextlib.contextmanager
def timer():
    t = []
    t.append(time.perf_counter())
    yield t
    t[-1] = time.perf_counter() - t[-1]

if __name__ == "__main__":
    hmx = pyhmx.HMx()
    camb.set_feedback_level(4)

    # Cosmological parameters for CAMB
    h = 0.7
    omc = 0.25
    omb = 0.048
    mnu = 0.0
    w = -1.0
    wa = 0.0
    ns = 0.97
    As = 2.1e-9

    z_max = 2.0
    n_z = 100

    k_max = 20.0

    halo_model_mode = pyhmx.constants.HMCode2016


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

    # Let CAMB do its thing
    with timer() as t_CAMB_lin:
        r = camb.get_results(p)

    # Get sigma_8, linear power spectrum, and Omega_m
    sigma8 = r.get_sigma8()[-1]
    # k_lin, z_lin, pofk_lin_camb = r.get_matter_power_spectrum(minkh=1e-3, maxkh=20.0, npoints=128)
    k_lin, z_lin, pofk_lin_camb = r.get_linear_matter_power_spectrum(nonlinear=False)

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
                  "As"   : 3.13}

    # Run HMCode
    with timer() as t_HMx:
        Pk_HMx_dmonly = hmx.run_HMCode(cosmology=cosmology,
                                       halo_model=halo_model,
                                       k=k_lin,
                                       z=z_lin,
                                       pk_lin=pofk_lin_camb,
                                       verbose=True)

    # Set up CAMB again (probably unnecessary)
    p = camb.CAMBparams(WantTransfer=True, 
                        WantCls=False, 
                        Want_CMB_lensing=False, 
                        DoLensing=False,
                        NonLinearModel=camb.nonlinear.Halofit(halofit_version="mead", HMCode_A_baryon=halo_model["As"], HMCode_eta_baryon=halo_model["eta0"]))
    p.set_cosmology(H0=h*100, omch2=omc*h**2, ombh2=omb*h**2, mnu=mnu)
    p.set_dark_energy(w=w, wa=wa)
    p.set_initial_power(camb.InitialPowerLaw(As=As, ns=ns))

    p.set_matter_power(redshifts=z_lin, kmax=max(k_lin), nonlinear=True)

    # Calculate CAMB outputs, now including HMCode
    with timer() as t_CAMB_nonlin:
        r = camb.get_results(p)

    Pk_nl_CAMB_interpolator = r.get_matter_power_interpolator()
    pofk_nonlin_camb = Pk_nl_CAMB_interpolator.P(z_lin, k_lin, grid=True)

    print("CAMB lin:     ", t_CAMB_lin[0])
    print("CAMB nonlin   ", t_CAMB_nonlin[0])
    print("CAMB diff:    ", t_CAMB_nonlin[0]- t_CAMB_lin[0])
    print("HMx:          ", t_HMx[0])