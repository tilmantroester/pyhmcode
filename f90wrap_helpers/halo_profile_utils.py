import numpy as np
import scipy.interpolate

import pyhmcode

try:
    import pyccl as ccl
except ImportError:
    print("halo_profile_utils provides utilities to interface pyhmcode with "
          "CCL. It appears that CCL is not installed, however.")
    raise


def ccl2hmcode_cosmo(ccl_cosmo, pofk_lin_z, pofk_lin_k_h, pofk_lin,
                     log10_T_heat=None):
    hmcode_cosmo = pyhmcode.Cosmology()

    hmcode_cosmo.om_m = ccl_cosmo["Omega_m"]
    hmcode_cosmo.om_b = ccl_cosmo["Omega_b"]
    hmcode_cosmo.om_v = ccl_cosmo["Omega_l"]
    hmcode_cosmo.h = ccl_cosmo["h"]
    hmcode_cosmo.ns = ccl_cosmo["n_s"]
    sigma8 = ccl_cosmo["sigma8"]
    if not np.isfinite(sigma8):
        sigma8 = ccl.sigma8(ccl_cosmo)
    hmcode_cosmo.sig8 = sigma8
    hmcode_cosmo.m_nu = ccl_cosmo["m_nu"].sum()

    if log10_T_heat is not None:
        hmcode_cosmo.theat = 10**log10_T_heat

    hmcode_cosmo.set_linear_power_spectrum(pofk_lin_k_h,
                                           pofk_lin_z,
                                           pofk_lin)

    return hmcode_cosmo


class HaloProfileInterpolated(ccl.halos.HaloProfile):
    def __init__(self, interpolator, is_logk=True, is_logM=True,
                 is_logp=True, norm=None):
        self.interpolator = interpolator
        self.is_logk = is_logk
        self.is_logM = is_logM
        self.is_logp = is_logp

        self.norm = norm or 1

    def _fourier(self, cosmo, k, M, a, mass_def):
        k_h = k/cosmo["h"]
        M_h = M*cosmo["h"]

        k_h = self._check_shape(k_h)
        M_h = self._check_shape(M_h)
        a = self._check_shape(a)

        k_h = np.log(k_h) if self.is_logk else k_h
        M_h = np.log(M_h) if self.is_logM else M_h

        positions = np.vstack([m.ravel() for m in np.meshgrid(a, k_h, M_h)]).T

        profile = self.interpolator(positions).reshape(len(k_h), len(M_h)).T

        if self.is_logp:
            profile = np.exp(profile)

        profile *= self.norm

        # Turn Mpc^3 h^-3 into Mpc^3
        profile *= cosmo["h"]**-3

        return profile.squeeze()

    def _check_shape(self, arr):
        if not isinstance(arr, np.ndarray):
            arr = np.array(arr)
        if arr.ndim < 1:
            arr = arr[None]
        return arr


class HMxProfileGenerator:
    def __init__(self, hmcode_cosmo, a_arr, k_arr, fields=None,
                 add_diffuse=False, verbose=False):
        self.hmcode_cosmo = hmcode_cosmo
        self.hmod = pyhmcode.Halomodel(
                        pyhmcode.HMx2020_matter_pressure_w_temp_scaling,
                        verbose=verbose)

        self.a_arr = a_arr.copy()
        self.k_arr = k_arr.copy()

        self.verbose = verbose

        if fields is None:
            self.fields = [pyhmcode.field_matter,
                           pyhmcode.field_electron_pressure]
        else:
            self.fields = fields

        self.add_diffuse = add_diffuse

        self._has_interpolator = False
        self._compute_tables()

    def _compute_tables(self):
        if self.verbose:
            print("Computing profile look-up table")

        profile_lut = np.empty((len(self.fields),
                                len(self.a_arr),
                                len(self.k_arr),
                                self.hmod.n), dtype=np.float64)

        for a_idx, a in enumerate(self.a_arr):
            if self.verbose:
                print("Computing a = ", a)
            pyhmcode.hmx.init_halomod(a, self.hmod,
                                      self.hmcode_cosmo, self.verbose)
            for k_idx, k in enumerate(self.k_arr):
                wk = np.empty((self.hmod.n, len(self.fields)), order="f")
                pyhmcode.hmx.init_windows(k, self.fields, len(self.fields),
                                          wk, self.hmod.n,
                                          self.hmod, self.hmcode_cosmo)
                if self.add_diffuse:
                    pyhmcode.hmx.add_smooth_component_to_windows(
                                        self.fields, len(self.fields),
                                        wk, self.hmod.n,
                                        self.hmod, self.hmcode_cosmo)
                profile_lut[:, a_idx, k_idx] = wk.T

        if pyhmcode.field_electron_pressure in self.fields:
            # Turn electron pressure into physical units
            pressure_idx = self.fields.index(pyhmcode.field_electron_pressure)
            profile_lut[pressure_idx] *= self.a_arr[:, None, None]**-3

        profile_lut[profile_lut <= 1e-30] = 1e-30

        profile_interpolator = [scipy.interpolate.RegularGridInterpolator(
                                            points=(self.a_arr,
                                                    np.log(self.k_arr),
                                                    np.log(self.hmod.m)),
                                            values=np.log(profile_lut[f_idx]),
                                            method="linear",
                                            bounds_error=True, fill_value=None)
                                for f_idx, f in enumerate(self.fields)]

        self.profiles = {f: HaloProfileInterpolated(
                                    profile_interpolator[f_idx],
                                    is_logk=True, is_logM=True, is_logp=True,
                                    norm=self.field_normalisation(f))
                         for f_idx, f in enumerate(self.fields)}

        self._has_interpolator = True

    def field_normalisation(self, field):
        if field == pyhmcode.field_cdm:
            return self.hmcode_cosmo.om_m/self.hmcode_cosmo.om_c
        else:
            return 1.0

    @property
    def matter_profile(self):
        if pyhmcode.field_matter in self.fields:
            return self.profiles[pyhmcode.field_matter]
        else:
            raise ValueError("Matter profiles not computed")

    @property
    def cdm_profile(self):
        if pyhmcode.field_cdm in self.fields:
            return self.profiles[pyhmcode.field_cdm]
        else:
            raise ValueError("CDM profiles not computed")

    @property
    def pressure_profile(self):
        if pyhmcode.field_electron_pressure in self.fields:
            return self.profiles[pyhmcode.field_electron_pressure]
        else:
            raise ValueError("Pressure profiles not computed")
