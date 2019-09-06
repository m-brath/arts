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
#include "doit.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "m_general.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
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
void newDOAngularGridsSet(  // WS Output:
    Index& doit_za_grid_size,
    Vector& scat_aa_grid,
    Vector& scat_za_grid,
    // Keywords:
    const Index& N_za_grid,
    const Index& N_aa_grid,
    const String& za_grid_opt_file,
    const Verbosity& verbosity) {
  // Azimuth angle grid (the same is used for the scattering integral and
  // for the radiative transfer.
  if (N_aa_grid > 1)
    nlinspace(scat_aa_grid, 0, 360, N_aa_grid);
  else if (N_aa_grid < 1) {
    ostringstream os;
    os << "N_aa_grid must be > 0 (even for 1D / DISORT cases).";
    throw runtime_error(os.str());
  } else {
    scat_aa_grid.resize(1);
    scat_aa_grid[0] = 0.;
  }

  // Zenith angle grid:
  // Number of zenith angle grid points (only for scattering integral):
  if (N_za_grid < 0) {
    ostringstream os;
    os << "N_za_grid must be >= 0.";
    throw runtime_error(os.str());
  } else
    doit_za_grid_size = N_za_grid;

  if (za_grid_opt_file == "")
    if (N_za_grid == 0)
      scat_za_grid.resize(0);
    else if (N_za_grid == 1) {
      ostringstream os;
      os << "N_za_grid must be >1 or =0 (the latter only allowed for RT4).";
      throw runtime_error(os.str());
    } else
      nlinspace(scat_za_grid, 0, 180, N_za_grid);
  else
    xml_read_from_file(za_grid_opt_file, scat_za_grid, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DoitCalc(Workspace& ws,
              Tensor7& doit_i_field,
              const Index& atmfields_checked,
              const Index& atmgeom_checked,
              const Index& cloudbox_checked,
              const Index& scat_data_checked,
              const Index& cloudbox_on,
              const Vector& f_grid,
              const Agenda& doit_mono_agenda,
              const Index& doit_is_initialized,
              const Verbosity& verbosity)

{
  CREATE_OUT2;

  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DOIT calculation will be skipped.\n";
    return;
    //throw runtime_error( "Cloudbox is off, no scattering calculations to be"
    //                     "performed." );
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

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;

  if (scat_data_checked != 1)
    throw runtime_error(
        "The scattering data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  chk_not_empty("doit_mono_agenda", doit_mono_agenda);

  // Frequency grid
  //
  if (f_grid.empty()) throw runtime_error("The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  // Check whether DoitInit was executed
  if (!doit_is_initialized) {
    ostringstream os;
    os << "Initialization method *DoitInit* has to be called before "
       << "*DoitCalc*";
    throw runtime_error(os.str());
  }

  //-------- end of checks ----------------------------------------

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws(ws);
  Agenda l_doit_mono_agenda(doit_mono_agenda);

  // OMP likes simple loop end conditions, so we make a local copy here:
  const Index nf = f_grid.nelem();

  if (nf) {
    String fail_msg;
    bool failed = false;

#pragma omp parallel for if (!arts_omp_in_parallel() && nf > 1) \
    firstprivate(l_ws, l_doit_mono_agenda)
    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) {
        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        continue;
      }

      try {
        ostringstream os;
        os << "Frequency: " << f_grid[f_index] / 1e9 << " GHz \n";
        out2 << os.str();

        Tensor6 doit_i_field_mono_local =
            doit_i_field(f_index, joker, joker, joker, joker, joker, joker);
        doit_mono_agendaExecute(
            l_ws, doit_i_field_mono_local, f_grid, f_index, l_doit_mono_agenda);
        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) =
            doit_i_field_mono_local;
      } catch (const std::exception& e) {
        doit_i_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        ostringstream os;
        os << "Error for f_index = " << f_index << " (" << f_grid[f_index]
           << " Hz)" << endl
           << e.what();
#pragma omp critical(DoitCalc_fail)
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


