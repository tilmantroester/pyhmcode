# pyhmcode

This library provides a convenient interface to Alexander Mead's library, which underlies the HMCode and HMx codes to predict non-linear power spectra. For details on the Fortran version, refer to https://github.com/alexander-mead/HMcode.
The `pyhmcode` interface uses the excellent [f90wrap](https://github.com/jameskermode/f90wrap) library to generate the python interface. This allows access to virtually all functionality in the library, although to speed up compilation, the interface is limited to commonly-used subroutines.
If you use `pyhmcode`, please cite [Tröster, Mead et al. 2021]() and the relevant papers describing the model:
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

## Citations
