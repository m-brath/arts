import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sensor

ws.frequency_grid = np.linspace(-50e6, 50e6, 1001) + 118750348044.712

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66-118e9-119e9"])
ws.ReadCatalogData()
ws.absorption_bandsSetZeeman(isot="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.options.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
alts = np.linspace(0.0, 1e5, 101)

ws.propagation_pathGeometric(pos=pos, los=los, max_step=alts[1] - alts[0])
ws.spectral_radiance_backgroundAgendasAtEndOfPath()
ws.propagation_path_atmospheric_pointFromPath()
ws.propagation_path_frequency_gridFromPath()
ws.propagation_path_propagation_matrixFromPath()
ws.propagation_path_transmission_matrixFromPath()
ws.propagation_path_transmission_matrix_cumulativeForward()
ws.propagation_path_spectral_radiance_sourceFromPropmat()
ws.propagation_path_spectral_radianceCalcEmission()
ws.background_transmittanceFromPathPropagationBack()
ws.spectral_radianceFromPathPropagation()


# %% Test calculations using single frequency approach
ws.spectral_radiance_operatorClearsky1D(altitude_grid=alts)
srad = ws.spectral_radiance_operator.geometric_planar(
    ws.frequency_grid, pos, [180 - los[0], 180 - los[1]]
)
assert np.allclose(srad, ws.spectral_radiance)
