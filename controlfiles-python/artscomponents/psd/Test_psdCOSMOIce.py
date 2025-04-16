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
            8.05178981e-26,
            1.06434989e-20,
            5.17585291e-16,
            9.25944489e-12,
            6.09387434e-08,
            1.47539245e-04,
            1.31409613e-01,
            4.30578394e01,
            5.19018513e03,
            2.30154242e05,
            3.75457371e06,
            2.25324183e07,
            4.97462667e07,
            4.04034732e07,
            1.20720889e07,
            1.32694116e06,
            5.36570050e04,
            7.98190779e02,
            4.36809930e00,
            8.79394647e-03,
            6.51299359e-06,
            1.77452859e-09,
            1.77865153e-13,
            6.55849601e-18,
            8.89658155e-23,
        ]
    ]
)

dpsd_data_dx_ref = np.array(
    [
        [
            [
                -5.99948439e-20,
                -7.31427593e-15,
                -3.25716442e-10,
                -5.29079222e-06,
                -3.12913241e-02,
                -6.72162600e01,
                -5.22584902e04,
                -1.46297785e07,
                -1.46292868e09,
                -5.15450031e10,
                -6.23456683e11,
                -2.43680572e12,
                -2.49928857e12,
                3.09700889e11,
                7.91580576e11,
                1.63846802e11,
                9.73246839e09,
                1.90998216e08,
                1.29817666e06,
                3.12273734e03,
                2.68991013e00,
                8.35647967e-04,
                9.40583983e-08,
                3.84803038e-12,
                5.73500732e-17,
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
