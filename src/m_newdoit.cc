/* Copyright (C) 2019

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. 
*/

/*!
  \file   m_newdoit.cc
  \author Manfred Brath manfred.brath@uni-hamburg.de
  \date   06/09/2019
  
  \brief  This file will contain the workingspace methods
          to calculate the radiative transfer
          inside the cloudbox using the newDoit method.

*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "m_general.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "newdoit.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "wsv_aux.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void NewDoitCalc(Workspace& ws,
                 // WS Output:
                 Tensor7& cloudbox_field,
                 Tensor7& cloudbox_field_clearsky,
                 Vector& za_grid,
                 Vector& aa_grid,
                 ArrayOfTensor3& gas_extinction_doit,
                 ArrayOfVector& p_grid_gas_extinction,
                 Tensor7& extinction_matrix_doit,
                 Tensor6& absorption_vector_doit,
                 ArrayOfTensor7& scattering_matrix_doit_array,
                 ArrayOfIndex& convergence_flag_doit,
                 ArrayOfIndex& iteration_counter_doit,


                 // WS Input
                 const Index& atmfields_checked,
                 const Index& atmgeom_checked,
                 const Index& scat_data_checked,
                 const Index& cloudbox_checked,
                 const Index& cloudbox_on,
                 const ArrayOfIndex& cloudbox_limits,
                 const Agenda& propmat_clearsky_agenda,
                 const Agenda& surface_rtprop_agenda,
                 const Agenda& iy_main_agenda,
                 const Agenda& ppath_step_agenda,
                 const Index& stokes_dim,
                 const Index& atmosphere_dim,
                 const Tensor4& pnd_field,
                 const Tensor3& t_field,
                 const EnergyLevelMap& nlte_field,
                 const Tensor3& z_field,
                 const Tensor4& vmr_field,
                 const Matrix& z_surface,
                 const Vector& p_grid,
                 const Vector& lat_grid,
                 const Vector& lon_grid,
                 const Vector& f_grid,
                 const Numeric& ppath_lmax,
                 const Numeric& ppath_lraytrace,
                 const ArrayOfArrayOfSingleScatteringData& scat_data,
                 const String& iy_unit,
                 const Vector& refellipsoid,

                 // Generic inputs
                 const Vector& epsilon,
                 const Index& max_num_iterations,
                 const Numeric& refinement_criterion,
                 const Numeric& tau_min,
                 const Index& N_decade,
                 const Index& accelerated,
                 const Index& ForwardCorrectionFlag,
                 const Index& t_interp_order,
                 const Index& N_za_grid,
                 const Index& N_aa_grid,
                 const Index& N_scat_za_grid,
                 const Index& N_scat_aa_grid,
                 const String& za_grid_type,
                 const Tensor7& cloudbox_field_apriori,
                 const Verbosity& verbosity) {
  CREATE_OUT1;
  CREATE_OUT2;

  if (!cloudbox_on) {
    out1 << "  Cloudbox is off, NewDoit calculation is skipped.\n";
    return;
  }

  //-------- Check input -------------------------------------------

  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  if (atmgeom_checked != 1)
    throw runtime_error(
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  if (cloudbox_checked != 1)
    throw runtime_error(
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  if (scat_data_checked != 1)
    throw runtime_error(
        "The scattering data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  if (!(atmosphere_dim == 1 || atmosphere_dim == 3))
    throw runtime_error("The atmospheric dimensionality must be 1 or 3.");

  // Frequency grid
  if (f_grid.empty()) throw runtime_error("The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  //-------- end of checks ----------------------------------------

  Vector scat_za_grid;
  Vector scat_aa_grid;

  SetAngularGrids(za_grid,
                  aa_grid,
                  scat_za_grid,
                  scat_aa_grid,
                  N_za_grid,
                  N_aa_grid,
                  N_scat_za_grid,
                  N_scat_aa_grid,
                  za_grid_type);

  //Check cloudbox_field_apriori
  if (!is_size(cloudbox_field_apriori, 0, 0, 0, 0, 0, 0, 0)) {
    chk_size("cloudbox_field_apriori",
             cloudbox_field_apriori,
             f_grid.nelem(),
             (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
             (cloudbox_limits[3] - cloudbox_limits[2]) + 1,
             (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
             N_za_grid,
             N_aa_grid,
             stokes_dim);

    cloudbox_field = cloudbox_field_apriori;
  } else {
    Initialize_cloudbox_field(cloudbox_field,
                            stokes_dim,
                            atmosphere_dim,
                            f_grid,
                            za_grid,
                            aa_grid,
                            cloudbox_limits);

    GetIncomingRadiation(ws,
                         //output
                         cloudbox_field,
                         //input
                         iy_main_agenda,
                         atmosphere_dim,
                         lat_grid,
                         lon_grid,
                         z_field,
                         nlte_field,
                         cloudbox_limits,
                         f_grid,
                         za_grid,
                         aa_grid,
                         verbosity);

    SetClearsky_cloudbox(cloudbox_field,
                         f_grid,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         cloudbox_limits,
                         atmosphere_dim,
                         verbosity);
  }

  for (Index l = 0; l < cloudbox_field.nlibraries(); l++)
    for (Index v = 0; v < cloudbox_field.nvitrines(); v++)
      for (Index s = 0; s < cloudbox_field.nshelves(); s++)
        for (Index b = 0; b < cloudbox_field.nbooks(); b++)
          for (Index p = 0; p < cloudbox_field.npages(); p++)
            for (Index r = 0; r < cloudbox_field.nrows(); r++)
              for (Index c = 0; c < cloudbox_field.ncols(); c++)
                if (std::isnan(cloudbox_field(l, v, s, b, p, r, c)))
                  throw std::runtime_error(
                      "*cloudbox_field_mono* contains at least one NaN value.\n"
                      "This indicates an improper initialization of *cloudbox_field*.");

  ostringstream os, os1, os2, os3;
  os << "limit fields to cloudbox \n";
  out1 << os.str();
  os.clear();

  Vector p_grid_cldbx;
  Vector lat_grid_cldbx;
  Vector lon_grid_cldbx;
  Tensor3 t_field_cldbx;
  Tensor3 z_field_cldbx;
  Tensor4 vmr_field_cldbx;

  cloudbox_field_clearsky=cloudbox_field;

  LimitInputGridsAndFieldsToCloudbox(p_grid_cldbx,
                                     lat_grid_cldbx,
                                     lon_grid_cldbx,
                                     t_field_cldbx,
                                     z_field_cldbx,
                                     vmr_field_cldbx,
                                     p_grid,
                                     lat_grid,
                                     lon_grid,
                                     t_field,
                                     z_field,
                                     vmr_field,
                                     cloudbox_limits,
                                     verbosity);



  os1 << "Now do loop... \n";
  out1 << os1.str();
  os1.clear();


  const Index nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;


  // Resize and initialize fields
  gas_extinction_doit.resize(nf);
  p_grid_gas_extinction.resize(nf);
  if (atmosphere_dim == 1) {

    extinction_matrix_doit.resize(
        nf, Np_cloud, 1, 1, za_grid.nelem(), stokes_dim, stokes_dim);
    absorption_vector_doit.resize(
        nf, Np_cloud, 1, 1, za_grid.nelem(), stokes_dim);

  } else {
    const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;

    extinction_matrix_doit.resize(nf,
                                  Np_cloud,
                                  Nlat_cloud,
                                  Nlon_cloud,
                                  za_grid.nelem(),
                                  stokes_dim,
                                  stokes_dim);
    absorption_vector_doit.resize(
        nf, Np_cloud, Nlat_cloud, Nlon_cloud, za_grid.nelem(), stokes_dim);
  }
  scattering_matrix_doit_array.resize(f_grid.nelem());
  convergence_flag_doit.resize(f_grid.nelem());
  iteration_counter_doit.resize(f_grid.nelem());

  extinction_matrix_doit = NAN;
  absorption_vector_doit = NAN;

  //now do loop over frequency to run DOIT for each frequency

  if (nf) {
    String fail_msg;
    bool failed = false;

    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) {
        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        gas_extinction_doit[f_index] = Tensor3();
        p_grid_gas_extinction[f_index] = Vector();
        extinction_matrix_doit(
            f_index, joker, joker, joker, joker, joker, joker) = 0;
        absorption_vector_doit(f_index, joker, joker, joker, joker, joker) =
            0;
        continue;
      }

      Tensor7 scattering_matrix;

      try {

        os2 << "Start with frequency: " << f_grid[f_index] / 1e9 << " GHz \n";
        out2 << os2.str();
        os2.clear();

        Tensor6 cloudbox_field_mono_local =
            cloudbox_field(f_index, joker, joker, joker, joker, joker, joker);

        Tensor3 gas_extinction_doit_local;
        Vector p_grid_gas_extinction_local;

        Tensor6 extinction_matrix_doit_local = extinction_matrix_doit(
            f_index, joker, joker, joker, joker, joker, joker);

        Tensor5 absorption_vector_doit_local =
            absorption_vector_doit(f_index, joker, joker, joker, joker, joker);

        Tensor7 scattering_matrix_doit_local;
        Index convergence_flag_local;
        Index iteration_counter_local;

        NewDoitMonoCalc(ws,
                        cloudbox_field_mono_local,
                        gas_extinction_doit_local,
                        extinction_matrix_doit_local,
                        absorption_vector_doit_local,
                        scattering_matrix_doit_local,
                        convergence_flag_local,
                        iteration_counter_local,
                        cloudbox_limits,
                        propmat_clearsky_agenda,
                        surface_rtprop_agenda,
                        ppath_step_agenda,
                        atmosphere_dim,
                        stokes_dim,
                        pnd_field,
                        t_field_cldbx,
                        z_field_cldbx,
                        vmr_field_cldbx,
                        z_surface,
                        p_grid_cldbx,
                        lat_grid_cldbx,
                        lon_grid_cldbx,
                        za_grid,
                        aa_grid,
                        scat_za_grid,
                        scat_aa_grid,
                        f_grid[f_index],
                        scat_data,
                        t_interp_order,
                        iy_unit,
                        refellipsoid,
                        epsilon,
                        max_num_iterations,
                        refinement_criterion,
                        tau_min,
                        N_decade,
                        accelerated,
                        ForwardCorrectionFlag,
                        ppath_lmax,
                        ppath_lraytrace,
                        verbosity);


        os3 << "Done with frequency: " << f_grid[f_index] / 1e9 << " GHz \n";
        out2 << os3.str();
        os3.clear();

        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) =
            cloudbox_field_mono_local;

        gas_extinction_doit[f_index] = gas_extinction_doit_local;

        p_grid_gas_extinction[f_index] = p_grid_gas_extinction_local;

        extinction_matrix_doit(f_index, joker, joker, joker, joker, joker, joker) =
            extinction_matrix_doit_local;

        absorption_vector_doit(f_index, joker, joker, joker, joker, joker) =
            absorption_vector_doit_local;

        scattering_matrix_doit_array[f_index]=scattering_matrix_doit_local;

        convergence_flag_doit[f_index]=convergence_flag_local;

        iteration_counter_doit[f_index]=iteration_counter_local;


      } catch (const std::exception& e) {
        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        gas_extinction_doit[f_index] = Tensor3();
        p_grid_gas_extinction[f_index] = Vector();
        extinction_matrix_doit(
            f_index, joker, joker, joker, joker, joker, joker) = NAN;
        absorption_vector_doit(
            f_index, joker, joker, joker, joker, joker) = NAN;


        os << "Error for f_index = " << f_index << " (" << f_grid[f_index]
           << " Hz)" << endl
           << e.what();
        {
          failed = true;
          fail_msg = os.str();
        }
        continue;
      }
    }

    if (failed) throw runtime_error(fail_msg);
  }
}
