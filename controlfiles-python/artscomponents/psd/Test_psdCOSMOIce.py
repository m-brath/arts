#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 11:11:23 2025

@author: Manfred Brath
"""
import numpy as np
import pyarts as pa


# Start workspace
ws = pa.Workspace()

# %% atmospheric condition

# Temperature
T = 250.0 # K

# Ice water content
IWC = 1e-5  # kg m^{-3}

# PSD grid in terms of maximum diameter
dmax = np.logspace(-7, -1, 25)

# %% the reference values

psd_data_ref = np.array(
    [
        [
            7.69880090e-30,
            1.01768899e-24,
            4.94894452e-20,
            8.85351261e-16,
            5.82672007e-12,
            1.41071153e-08,
            1.25648641e-05,
            4.11701920e-03,
            4.96264842e-01,
            2.20064324e01,
            3.58997393e02,
            2.15446015e03,
            4.75653999e03,
            3.86321927e03,
            1.15428508e03,
            1.26876832e02,
            5.13046922e00,
            7.63198249e-02,
            4.17660267e-04,
            8.40842156e-07,
            6.22746521e-10,
            1.69673360e-13,
            1.70067579e-17,
            6.27097281e-22,
            8.50655713e-27,
        ]
    ]
)

dpsd_data_dx_ref = np.array(
    [
        [
            [
                -5.73646815e-24,
                -6.99361949e-19,
                -3.11437096e-14,
                -5.05884491e-10,
                -2.99195185e-06,
                -6.42695122e-03,
                -4.99674881e00,
                -1.39884118e03,
                -1.39879417e05,
                -4.92852801e06,
                -5.96124462e07,
                -2.32997663e08,
                -2.38972024e08,
                2.96123662e07,
                7.56877837e07,
                1.56663790e07,
                9.30579886e05,
                1.82624891e04,
                1.24126484e02,
                2.98583714e-01,
                2.57198499e-04,
                7.99013321e-08,
                8.99348963e-12,
                3.67933347e-16,
                5.48358570e-21,
            ]
        ]
    ]
)


# %%  psd and dspd

ws.psd_size_grid = dmax
ws.pnd_agenda_input_t = pa.arts.Vector([T])
ws.pnd_agenda_input = pa.arts.Matrix(np.array([[IWC]]))
ws.pnd_agenda_input_names = pa.arts.ArrayOfString(["IWC-mass_density"])
ws.dpnd_data_dx_names = pa.arts.ArrayOfString(["IWC-mass_density"])
ws.psdCOSMOIce()

# Comapre the results with the reference values
ws.CompareRelative(ws.psd_data, psd_data_ref, 1e-6)
ws.CompareRelative(ws.dpsd_data_dx, dpsd_data_dx_ref, 1e-6)
