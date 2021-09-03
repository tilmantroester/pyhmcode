import os
import ctypes as ct
import numpy as np

def _array_ctype(ndim, dtype=np.float64, flags="F"):
    return [ct.POINTER(ct.c_int)]*ndim + [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)]

def _array_arg(a):
    arr = a
    return (*(ct.c_int(s) for s in arr.shape), arr)

def _load_lib(libname, path=None):
    if path is None:
        path = os.path.dirname(__file__)
    libpath = os.path.abspath(os.path.join(path, libname))
    return ct.CDLL(libpath)

_hmx_lib = _load_lib("libhmx_wrapper.so")

class HMxConstants:
    def __init__(self, lib):
        for c in ["HMCode2016", "HMCode2016_CAMB", "HMCode2020", 
                  "HMx2020_matter_with_temperature_scaling", "HMx2020_matter_pressure_with_temperature_scaling",
                  "field_dmonly", "field_matter", "field_cdm", "field_gas", "field_stars", "field_electron_pressure"]:
            setattr(self, c, ct.c_int.in_dll(lib, f"constant_{c.lower()}").value)

constants = HMxConstants(_hmx_lib)


class HMx:
    module_name = "HMx_wrapper"

    def __init__(self):
        self.lib = _hmx_lib

    def run_HMCode(self, cosmology=None, halo_model=None, 
                   k=None, z=None, 
                   pk_lin=None,
                   verbose=False):
        cosmology = cosmology or {}
        halo_model = halo_model or {}

        pofk = self._run_hmx(cosmology.get("Omega_m"), cosmology.get("Omega_b"), cosmology.get("Omega_v"), 
                             cosmology.get("h"), cosmology.get("n_s"), cosmology.get("sigma_8"), 
                             cosmology.get("m_nu", 0.06), cosmology.get("w", -1.0), cosmology.get("w_a", 0.0), 
                             halo_model_mode=constants.HMCode2016, Theat=10**7.8, 
                             eta0=halo_model.get("eta0", 0.603), As=halo_model.get("A", 3.13),
                             fields=np.array([constants.field_dmonly]),
                             k=k, z=z,
                             pk_lin=pk_lin,
                             verbose=verbose)
        return pofk[0,0]

    def run_HMx(self, cosmology=None, halo_model=None, fields=None,
                mode=constants.HMx2020_matter_with_temperature_scaling,
                k=None, z=None, pk_lin=None,
                verbose=False):
        cosmology = cosmology or {}
        halo_model = halo_model or {}
        fields = fields or [constants.field_matter]
        fields = np.array(fields)

        pofk = self._run_hmx(cosmology.get("Omega_m"), cosmology.get("Omega_b"), cosmology.get("Omega_v"), 
                             cosmology.get("h"), cosmology.get("n_s"), cosmology.get("sigma_8"), 
                             cosmology.get("m_nu", 0.06), cosmology.get("w", -1.0), cosmology.get("w_a", 0.0), 
                             halo_model_mode=mode, Theat=halo_model.get("Theat", 10**7.8),
                             eta0=0.603, As=3.13,
                             fields=fields,
                             k=k, z=z,
                             pk_lin=pk_lin,
                             verbose=verbose)
        return pofk

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def _run_hmx(self, omm, omb, omv, h, ns, sigma8, mnu, w, wa,
                       halo_model_mode, Theat, eta0, As,
                       fields=None,
                       k=None, z=None,
                       pk_lin=None,
                       verbose=True):
        f = self.get_function("run_HMx")
        f.restype = ct.c_int
        f.argtypes = [ct.POINTER(ct.c_double),     # omm
                      ct.POINTER(ct.c_double),     # omb
                      ct.POINTER(ct.c_double),     # omv
                      ct.POINTER(ct.c_double),     # mnu
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # ns
                      ct.POINTER(ct.c_double),     # sigma8
                      ct.POINTER(ct.c_double),     # w
                      ct.POINTER(ct.c_double),     # wa
                      ct.POINTER(ct.c_int),        # halo_model_mode
                      ct.POINTER(ct.c_double),     # Theat
                      ct.POINTER(ct.c_double),     # eta0
                      ct.POINTER(ct.c_double),     # As
                      *_array_ctype(ndim=1, dtype=np.int32), # fields
                      *_array_ctype(ndim=1, dtype=np.float64), # k
                      *_array_ctype(ndim=1, dtype=np.float64), # a
                      *_array_ctype(ndim=2, dtype=np.float64), # Pk_lin
                      *_array_ctype(ndim=4, dtype=np.float64), # Pk_hmx
                      ct.POINTER(ct.c_bool),        # verbose
                      ]      

        Pk_hmx = np.zeros((len(fields), len(fields), len(k), len(z)), dtype=np.float64, order="F")

        if k is None or z is None or pk_lin is None:
            raise ValueError("k, z, and pk_lin need to be specified.")

        if (len(z), len(k)) != pk_lin.shape:
            raise ValueError("Shape of pk_lin does not match z and k arrays.")

        if len(z) > 1 and z[0] > z[1]:
            raise ValueError("Redshift needs to be increasing.")

        a = 1/(1+np.array(z))

        status = f(ct.c_double(omm), ct.c_double(omb), ct.c_double(omv), ct.c_double(mnu), 
          ct.c_double(h), ct.c_double(ns), ct.c_double(sigma8), ct.c_double(w), ct.c_double(wa),
          ct.c_int(halo_model_mode), ct.c_double(Theat), ct.c_double(eta0), ct.c_double(As),
          *_array_arg(np.ascontiguousarray(fields, dtype=np.int32)),
          *_array_arg(np.ascontiguousarray(k, dtype=np.float64)),
          *_array_arg(np.ascontiguousarray(a[::-1], dtype=np.float64)),      # Reverse order for HMx
          *_array_arg(np.asfortranarray(pk_lin[::-1].T, dtype=np.float64)),  # Reverse order and transpose to (k, z) for HMx
          *_array_arg(Pk_hmx),
          ct.c_bool(verbose)
          )
        if status != 0:
            raise RuntimeError("HMx failed.")
        # Restore CAMB order
        return np.swapaxes(Pk_hmx[...,::-1], 2, 3)

