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
                 Tensor7& doit_i_field,
                 Vector& za_grid,
                 Vector& aa_grid,
                 Vector& scat_za_grid,
                 Vector& scat_aa_grid,
                 Tensor4& gas_extinct_doit,

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
                 const Tensor4& nlte_field,
                 const Tensor3& z_field,
                 const Tensor4& vmr_field,
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
                 const Index& max_lvl_optimize,
                 const Numeric& tau_scat_max,
                 const Numeric& sgl_alb_max,
                 const Index& N_za_grid,
                 const Index& N_aa_grid,
                 const Index& N_scat_za_grid,
                 const Index& N_scat_aa_grid,
                 const String& za_grid_type,
                 const Tensor7& doit_i_field_apriori,
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

  SetAngularGrids(za_grid,
                  aa_grid,
                  scat_za_grid,
                  scat_aa_grid,
                  N_za_grid,
                  N_aa_grid,
                  N_scat_za_grid,
                  N_scat_aa_grid,
                  za_grid_type);

  //Check doit_i_field_apriori
  if (!is_size(doit_i_field_apriori, 0, 0, 0, 0, 0, 0, 0)) {
    chk_size("doit_i_field_apriori",
             doit_i_field_apriori,
             f_grid.nelem(),
             (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
             (cloudbox_limits[3] - cloudbox_limits[2]) + 1,
             (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
             N_za_grid,
             N_aa_grid,
             stokes_dim);

    doit_i_field = doit_i_field_apriori;
  } else {
    Initialize_doit_i_field(doit_i_field,
                            stokes_dim,
                            atmosphere_dim,
                            f_grid,
                            za_grid,
                            aa_grid,
                            cloudbox_limits);

    GetIncomingRadiation(ws,
                         //output
                         doit_i_field,
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

    SetClearsky_doit_i_field(doit_i_field,
                             f_grid,
                             p_grid,
                             lat_grid,
                             lon_grid,
                             cloudbox_limits,
                             atmosphere_dim,
                             verbosity);
  }

  ostringstream os;
  os << "limit fields to cloudbox \n";
  out1 << os.str();

  Vector p_grid_cldbx;
  Vector lat_grid_cldbx;
  Vector lon_grid_cldbx;
  Tensor3 t_field_cldbx;
  Tensor3 z_field_cldbx;
  Tensor4 vmr_field_cldbx;

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



  os << "Now do loop... \n";
  out1 << os.str();


  const Index nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  // Resize and initialize gas absorption field inside cloudbox
  if (atmosphere_dim == 1) {
    gas_extinct_doit.resize(nf, Np_cloud, 1, 1);
  } else {
    const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    gas_extinct_doit.resize(nf, Np_cloud, Nlat_cloud, Nlon_cloud);
  }
  gas_extinct_doit = NAN;

  //now do loop over frequency to run DOIT for each frequency

  if (nf) {
    String fail_msg;
    bool failed = false;

    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) {
        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        gas_extinct_doit(f_index, joker, joker, joker) = NAN;
        continue;
      }

      try {
        ostringstream os;
        os << "Start with frequency: " << f_grid[f_index] / 1e9 << " GHz \n";
        out2 << os.str();

        Tensor6 doit_i_field_mono_local =
            doit_i_field(f_index, joker, joker, joker, joker, joker, joker);

        Tensor3 gas_extinct_doit_local =
            gas_extinct_doit(f_index, joker, joker, joker);

        NewDoitMonoCalc(ws,
                        doit_i_field_mono_local,
                        gas_extinct_doit_local,
                        cloudbox_limits,
                        propmat_clearsky_agenda,
                        surface_rtprop_agenda,
                        atmosphere_dim,
                        stokes_dim,
                        pnd_field,
                        t_field_cldbx,
                        z_field_cldbx,
                        vmr_field_cldbx,
                        p_grid_cldbx,
                        lat_grid_cldbx,
                        lon_grid_cldbx,
                        za_grid,
                        aa_grid,
                        scat_za_grid,
                        scat_aa_grid,
                        f_grid[f_index],
                        scat_data,
                        iy_unit,
                        refellipsoid,
                        epsilon,
                        max_num_iterations,
                        max_lvl_optimize,
                        tau_scat_max,
                        sgl_alb_max,
                        verbosity);

        os << "Done with frequency: " << f_grid[f_index] / 1e9 << " GHz \n";
        out2 << os.str();

        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) =
            doit_i_field_mono_local;

        gas_extinct_doit(f_index, joker, joker, joker) = gas_extinct_doit_local;

      } catch (const std::exception& e) {
        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        gas_extinct_doit(f_index, joker, joker, joker) = -9999.;
        ostringstream os;
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
