# pyhmcode

This library provides a convenient interface to Alexander Mead's library, which underlies the HMCode and HMx codes to predict non-linear power spectra.
For details on the Fortran version, refer to https://github.com/alexander-mead/HMcode. 
If you use this code, please cite the relevant papers describing the model
- HMCode2015: Mead et al. (2015; https://arxiv.org/abs/1505.07833)
- HMCode2016: Mead et al. (2016; https://arxiv.org/abs/1602.02154)
- HMCode2020: Mead et al. (2021; https://arxiv.org/abs/2009.01858)
- HMx: Mead, Tröster et al. (2020; https://arxiv.org/abs/2005.00009)

# Installation
```
git clone --recursive https://github.com/alexander-mead/HMcode
```

```
pip install pyhmcode
```

# Alternative interface and CosmoSIS support
There is an alternative interface to compute powerspectra in the `powerspectrum_interface` subdirectory. 
This interface uses `ctypes` and a thin Fortran wrapper to access the main functionality of HMCode and HMx, namely predicting the non-linear power spectra.
Installation proceeds by 
```
cd powerspectrum_interface
pip install .
```
This installs the `pyhmx` python module, which in turn is interfaced with CosmoSIS using the `cosmosis_interface.py` module.

