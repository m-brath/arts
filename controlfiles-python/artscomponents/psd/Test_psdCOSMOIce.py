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
            2.97959835e-26,
            3.93867107e-21,
            1.91534592e-16,
            3.42649613e-12,
            2.25506357e-08,
            5.45975119e-05,
            4.86286744e-02,
            1.59337327e01,
            1.92064962e03,
            8.51695354e04,
            1.38939564e06,
            8.33821522e06,
            1.84088131e07,
            1.49514735e07,
            4.46732676e06,
            4.91040103e05,
            1.98559981e04,
            2.95373821e02,
            1.61643333e00,
            3.25423651e-03,
            2.41016040e-06,
            6.56671695e-10,
            6.58197406e-14,
            2.42699875e-18,
            3.29221704e-23,
        ]
    ]
)

dpsd_data_dx_ref = np.array(
    [
        [
            [
                -1.92217434e-20,
                -2.31281115e-15,
                -1.01379266e-10,
                -1.61523008e-06,
                -9.32442118e-03,
                -1.94139057e01,
                -1.44756044e04,
                -3.82043717e06,
                -3.49297887e08,
                -1.05574895e10,
                -9.17731777e10,
                -6.79285812e10,
                9.16009170e11,
                1.60975345e12,
                7.39660364e11,
                1.09736202e11,
                5.58714026e09,
                1.00217067e08,
                6.42039011e05,
                1.48100561e03,
                1.23642848e00,
                3.74902179e-04,
                4.13886760e-08,
                1.66667954e-12,
                2.45148504e-17,
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
