#include "workspace_agendas.h"

std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas_creator() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["propagation_matrix_agenda"] = {
      .desc =
          R"--(Compute the propagation matrix, the non-LTE source vector, and their derivatives
)--",
      .output       = {"propagation_matrix",
                       "propagation_matrix_source_vector_nonlte",
                       "propagation_matrix_jacobian",
                       "propagation_matrix_source_vector_nonlte_jacobian"},
      .input        = {"frequency_grid",
                       "frequency_grid_wind_shift_jacobian",
                       "jacobian_targets",
                       "select_species",
                       "ray_path_point",
                       "atmospheric_point"},
      .enum_options = {"Empty"},
  };

  wsa_data["propagation_matrix_scattering_spectral_agenda"] = {
      .desc =
          R"--(Get the scattering propagation matrix, the scattering absorption vector, and the scattering spectral phase matrix
)--",
      .output = {"propagation_matrix_scattering",
                 "absorption_vector_scattering",
                 "phase_matrix_scattering_spectral"},
      .input  = {"frequency_grid", "atmospheric_point", "legendre_degree"},
      .enum_options = {"FromSpeciesTRO"},
      .enum_default = "FromSpeciesTRO",
  };

  wsa_data["propagation_matrix_scattering_agenda"] = {
      .desc =
          R"--(Compute the propagation matrix, the non-LTE source vector, and their derivatives
)--",
      .output       = {"propagation_matrix_scattering"},
      .input        = {"frequency_grid", "atmospheric_point"},
      .enum_options = {"AirSimple"},
      .enum_default = "AirSimple",
  };

  wsa_data["ray_path_observer_agenda"] = {
      .desc         = R"--(Get the propagation path as it is obeserved.

The intent of this agenda is to provide a propagation path as seen from the observer
position and line of sight.
)--",
      .output       = {"ray_path"},
      .input        = {"spectral_radiance_observer_position",
                       "spectral_radiance_observer_line_of_sight"},
  };

  wsa_data["spectral_radiance_observer_agenda"] = {
      .desc =
          R"--(Spectral radiance as seen from the input position and environment

The intent of this agenda is to provide a spectral radiance as seen from the observer
position and line of sight.

It also outputs the *ray_path* as seen from the observer position and line of sight.
This is useful in-case a call to the destructive *spectral_radianceApplyUnitFromSpectralRadiance*
is warranted

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
- *ray_path* : (Unknown)
)--",
      .output = {"spectral_radiance", "spectral_radiance_jacobian", "ray_path"},
      .input  = {"frequency_grid",
                 "jacobian_targets",
                 "spectral_radiance_observer_position",
                 "spectral_radiance_observer_line_of_sight",
                 "atmospheric_field",
                 "surface_field"},
      .enum_options = {"Emission"},
      .enum_default = "Emission",
  };

  wsa_data["spectral_radiance_space_agenda"] = {
      .desc         = R"--(Spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at space.

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
)--",
      .output       = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input        = {"frequency_grid", "jacobian_targets", "ray_path_point"},
      .enum_options = {"UniformCosmicBackground",
                       "SunOrCosmicBackground",
                       "Transmission"},
      .enum_default = "UniformCosmicBackground",
  };

  wsa_data["spectral_radiance_surface_agenda"] = {
      .desc         = R"--(Spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at the surface.

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
)--",
      .output       = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input        = {"frequency_grid",
                       "jacobian_targets",
                       "ray_path_point",
                       "surface_field"},
      .enum_options = {"Blackbody", "Transmission"},
      .enum_default = "Blackbody",
  };

  wsa_data["inversion_iterate_agenda"] = {
      .desc   = R"--(Work in progress ...

The WSV *measurement_jacobian* is both in- and output. As input variable, *measurement_jacobian*
is assumed to be valid for the previous iteration. For the first iteration
the input *measurement_jacobian* shall be set to have size zero, to flag that there
is not yet any calculated Jacobian.
)--",
      .output = {"measurement_vector_fitted", "measurement_jacobian"},
      .input  = {"model_state_vector",
                 "inversion_iterate_agenda_do_jacobian",
                 "inversion_iterate_agenda_counter"},
  };

  wsa_data["disort_settings_agenda"] = {
      .desc         = R"--(An agenda for setting up Disort.

The only intent of this Agenda is to simplify the setup of Disort for different
scenarios.  The output of this Agenda is just that setting.
)--",
      .output       = {"disort_settings"},
      .input        = {"frequency_grid",
                       "ray_path",
                       "disort_quadrature_dimension",
                       "disort_fourier_mode_dimension",
                       "disort_legendre_polynomial_dimension"},
  };

  return wsa_data;
}

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>&
internal_workspace_agendas() {
  static const auto out = internal_workspace_agendas_creator();
  return out;
}
