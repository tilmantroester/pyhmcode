import os
import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import pyhmx

mode_mapping = {"HMCode2016"              : pyhmx.constants.HMCode2016,
                "HMx2020_matter"          : pyhmx.constants.HMx2020_matter_with_temperature_scaling,
                "HMx2020_matter_pressure" : pyhmx.constants.HMx2020_matter_pressure_with_temperature_scaling
               }

field_mapping = {"matter"              : pyhmx.constants.field_matter,
                 "cdm"                 : pyhmx.constants.field_cdm,
                 "stars"               : pyhmx.constants.field_stars,
                 "gas"                 : pyhmx.constants.field_gas,
                 "pressure"            : pyhmx.constants.field_electron_pressure,}

def setup(options):
    config = {}

    config["module"] = pyhmx.HMx()

    mode = options.get_string(option_section, "mode", default="HMCode2016")
    try:
        mode = mode_mapping[mode]
    except KeyError:
        raise ValueError(f"mode {mode} is not supported. Choose one of {list(mode_mapping.keys())}.")
    config["mode"] = mode

    field_strings = options.get_string(option_section, "fields", default="")
    if mode == pyhmx.constants.HMCode2016:
        fields = [pyhmx.constants.field_dmonly]
    else:
        field_strings = [s.lower().strip() for s in field_strings.split(" ")]
        fields = []
        for f in field_strings:
            try:
                fields.append(field_mapping[f])
            except KeyError:
                raise ValueError(f"Field {f} not supported. Supported fields are {list(field_mapping.keys())}.")
        
        config["field_names"] = field_strings

    config["fields"] = fields
    
    config["verbose"] = options.get_int(option_section, "verbose", 1)
    config["log10_Theat"] = options.get_bool(option_section, "log10_Theat", default=True)

    config["nonlinear_matter_matter_power_output_section"] = options.get_string(option_section, "nonlinear_matter_matter_power_output_section", names.matter_power_nl)
    config["linear_matter_power_input_section"] = options.get_string(option_section, "linear_matter_power_input_section", names.matter_power_lin)

    return config

def execute(block, config):
    cosmology = {"Omega_m"  : block[names.cosmological_parameters, "omega_m"],
                 "Omega_b"  : block[names.cosmological_parameters, "omega_b"],
                 "Omega_v"  : block[names.cosmological_parameters, "omega_lambda"],
                 "h"        : block[names.cosmological_parameters, "h0"],
                 "n_s"      : block[names.cosmological_parameters, "n_s"],
                 "sigma_8"  : block[names.cosmological_parameters, "sigma_8"],
                 "mnu"      : block[names.cosmological_parameters, "mnu"],
                 "w"        : block.get_double(names.cosmological_parameters, "w", -1.0),
                 "wa"       : block.get_double(names.cosmological_parameters, "wa", 0.0),
                 }

    pk_lin = block[names.matter_power_lin, "p_k"]
    k_h = block[names.matter_power_lin, "k_h"]
    z = block[names.matter_power_lin, "z"]

    if config["mode"] == pyhmx.constants.HMCode2016:
        # Compute matter power spectrum with HMCode, parametrised by eta0 and A
        halo_model = {"eta0" : block[names.halo_model_parameters, "eta0"],
                      "A"    : block[names.halo_model_parameters, "A"]}
        
        Pk_HMx_dmonly = config["module"].run_HMCode(cosmology=cosmology,
                                                    halo_model=halo_model,
                                                    k=k_h,
                                                    z=z,
                                                    pk_lin=pk_lin,
                                                    verbose=config["verbose"])

        block.put_grid(config["nonlinear_matter_matter_power_output_section"], "z", z, "k_h", k_h, "p_k", Pk_HMx_dmonly)
    else:
        # Use HMx to compute the power spectra
        if config["log10_Theat"]:
            halo_model = {"Theat" : 10**block[names.halo_model_parameters, "log10_Theat"]}
        else:
            halo_model = {"Theat" : block[names.halo_model_parameters, "Theat"]}

        Pk_HMx = config["module"].run_HMx(cosmology=cosmology, halo_model=halo_model,
                                          fields=config["fields"],
                                          mode=config["mode"],
                                          k=k_h,
                                          z=z,
                                          pk_lin=pk_lin,
                                          verbose=config["verbose"])

        for i, field_name_i in enumerate(config["field_names"]):
            for j, field_name_j in enumerate(config["field_names"][i:]):
                block.put_grid(f"{field_name_i}_{field_name_j}_power_spectrum", "z", z, "k_h", k_h, "p_k", Pk_HMx[i,j+i])

                if field_name_i == "matter" and field_name_j == "matter":
                    block.put_grid(config["nonlinear_matter_matter_power_output_section"], "z", z, "k_h", k_h, "p_k", Pk_HMx[i,j+i])

    return 0

def cleanup(config):
    pass
