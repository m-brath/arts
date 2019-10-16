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
  \file   newdoit.cc
  \author Manfred Brath manfred.brath@uni-hamburg.de
  \date   06/09/2019
  
  \brief  This file will contain the internal functions
          to calculate the radiative transfer
          inside the cloudbox using the newDoit method.
  
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "newdoit.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "propagationmatrix.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;


void InitializeDoitFields(  //Output
    Tensor7& doit_i_field,
    // WS Input
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const Vector& rt_za_grid,
    const Vector& rt_aa_grid,
    const ArrayOfIndex& cloudbox_limits) {
  //------------- end of checks ---------------------------------------

  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Nrza = rt_za_grid.nelem();
  const Index Ns = stokes_dim;

  // Resize and initialize radiation field in the cloudbox
  if (atmosphere_dim == 1) {
    doit_i_field.resize(Nf, Np_cloud, 1, 1, Nrza, 1, Ns);
  } else if (atmosphere_dim == 3) {
    const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    const Index Nraa = rt_aa_grid.nelem();

    doit_i_field.resize(Nf, Np_cloud, Nlat_cloud, Nlon_cloud, Nrza, Nraa, Ns);

  }

  doit_i_field = NAN;
}

void New_DOITAngularGridsSet(Vector &za_grid,
                          Vector &aa_grid,
                          Vector &scat_za_grid,
                          Vector &scat_aa_grid,
    // Keywords:
                          const Index &N_za_grid,
                          const Index &N_aa_grid,
                          const Index &N_scat_za_grid,
                          const Index &N_scat_aa_grid,
                          const String &za_grid_type) {
  // Azimuth angle grid
  if (N_aa_grid > 1)
    nlinspace(aa_grid, 0, 360, N_aa_grid);
  else if (N_aa_grid < 1) {
    ostringstream os;
    os << "N_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    aa_grid.resize(1);
    aa_grid[0] = 0.;
  }

  // Azimuth angle grid (scattering)
  if (N_scat_aa_grid > 1)
    nlinspace(scat_aa_grid, 0, 360, N_scat_aa_grid);
  else if (N_scat_aa_grid < 1) {
    ostringstream os;
    os << "N_scat_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    scat_aa_grid.resize(1);
    scat_aa_grid[0] = 0.;
  }

  // Make scattering zenith angle
  nlinspace(scat_za_grid, 0, 180, N_scat_za_grid);

  //calculate zenith angle grid
  za_grid.resize(N_za_grid);
  za_grid = 0.;
  if (za_grid_type == "linear") {
    nlinspace(za_grid, 0, 180, N_za_grid);


  } else if (za_grid_type == "linear_mu") {
    Vector x;
    nlinspace(x, -1, 1, N_za_grid);

    //allocate
    Vector za_grid_temp;
    za_grid_temp.resize(x.nelem());

    for (Index i = 0; i < N_za_grid; i++) {
      za_grid_temp[i] = acos(x[i]) / DEG2RAD;
    }

    //#sort weights and theta in increasing direction of scat_za_grid
    za_grid = za_grid_temp[Range(x.nelem() - 1, x.nelem(), -1)];


  } else {
    ostringstream os;
    os << "The selected grid type is not implemented";
    throw std::runtime_error(os.str());
  }

  //be sure that the first and the last angle are within the closed interval
  //between 0 and 180 deg, because ARTS is picky if the angles are due to numerics
  // slightly below and above,respectively.
  if (za_grid[0] < 0) {
    za_grid[0] = 0.;
  }

  if (scat_za_grid[0] < 0) {
    scat_za_grid[0] = 0.;
  }

  if (za_grid[za_grid.nelem() - 1] > 180) {
    za_grid[za_grid.nelem() - 1] = 180.;
  }

  if (scat_za_grid[scat_za_grid.nelem() - 1] > 180) {
    scat_za_grid[scat_za_grid.nelem() - 1] = 180.;
  }

}


void New_doit_i_fieldSetClearsky(Tensor7& doit_i_field,
                             const Vector& f_grid,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Vector& za_grid,
                             const Vector& aa_grid,
                             const ArrayOfIndex& cloudbox_limits,
                             const Index& atmosphere_dim,
                             const Index& stokes_dim,
                             const Verbosity& verbosity) {
  CREATE_OUT2;

  out2<< "  Interpolate boundary clearsky field to obtain the initial field.\n";


  InitializeDoitFields(doit_i_field,
      stokes_dim,
      atmosphere_dim,
      f_grid,
      za_grid,
      aa_grid,
      cloudbox_limits);



  // Initial field only needs to be calculated from clearsky field for the
  // first frequency. For the next frequencies the solution field from the
  // previous frequencies is used.
  if (atmosphere_dim == 1) {
    const Index nf = f_grid.nelem() ? doit_i_field.nlibraries() : 1;

    for (Index f_index = 0; f_index < nf; f_index++) {
      Index N_za = doit_i_field.npages();
      Index N_aa = doit_i_field.nrows();
      Index N_i = doit_i_field.ncols();

      //1. interpolation - pressure grid

      /*the old grid is having only two elements, corresponding to the
             cloudbox_limits and the new grid have elements corresponding to
             all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1] - cloudbox_limits[0]) + 1);

      p2gridpos(p_gp,
                p_grid[Range(cloudbox_limits[0],
                             2,
                             (cloudbox_limits[1] - cloudbox_limits[0]))],
                p_grid[Range(cloudbox_limits[0],
                             (cloudbox_limits[1] - cloudbox_limits[0]) + 1)]);

      Matrix itw((cloudbox_limits[1] - cloudbox_limits[0]) + 1, 2);
      interpweights(itw, p_gp);

      Tensor6 scat_i_p(2, 1, 1, N_za, 1, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          doit_i_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          doit_i_field(f_index,
                       doit_i_field.nvitrines() - 1,
                       joker,
                       joker,
                       joker,
                       joker,
                       joker);

      for (Index za_index = 0; za_index < N_za; ++za_index) {
        for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
          for (Index i = 0; i < N_i; ++i) {
            VectorView target_field = doit_i_field(
                f_index, Range(joker), 0, 0, za_index, aa_index, i);

            ConstVectorView source_field =
                scat_i_p(Range(joker), 0, 0, za_index, aa_index, i);

            interp(target_field, itw, source_field, p_gp);
          }
        }
      }
    }
  } else if (atmosphere_dim == 3) {
    if (f_grid.nelem() == 0)
      throw runtime_error(
          "Error in doit_i_fieldSetClearsky: For 3D "
          "all_frequencies option is not implemented \n");

    for (Index f_index = 0; f_index < doit_i_field.nvitrines(); f_index++) {
      Index N_p = doit_i_field.nvitrines();
      Index N_lat = doit_i_field.nshelves();
      Index N_lon = doit_i_field.nbooks();
      Index N_za = doit_i_field.npages();
      Index N_aa = doit_i_field.nrows();
      Index N_i = doit_i_field.ncols();

      Tensor6 scat_i_p(2, N_lat, N_lon, N_za, N_aa, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          doit_i_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          doit_i_field(f_index, N_p - 1, joker, joker, joker, joker, joker);

      Tensor6 scat_i_lat(N_p, 2, N_lon, N_za, N_aa, N_i);
      scat_i_lat(joker, 0, joker, joker, joker, joker) =
          doit_i_field(f_index, joker, 0, joker, joker, joker, joker);
      scat_i_lat(joker, 1, joker, joker, joker, joker) =
          doit_i_field(f_index, joker, N_lat - 1, joker, joker, joker, joker);

      Tensor6 scat_i_lon(N_p, N_lat, 2, N_za, N_aa, N_i);
      scat_i_lon(joker, joker, 0, joker, joker, joker) =
          doit_i_field(f_index, joker, joker, 0, joker, joker, joker);
      scat_i_lon(joker, joker, 1, joker, joker, joker) =
          doit_i_field(f_index, joker, joker, N_lon - 1, joker, joker, joker);

      //1. interpolation - pressure grid, latitude grid and longitude grid

      ArrayOfGridPos p_gp((cloudbox_limits[1] - cloudbox_limits[0]) + 1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3] - cloudbox_limits[2]) + 1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5] - cloudbox_limits[4]) + 1);

      /*the old grid is having only two elements, corresponding to the
             cloudbox_limits and the new grid have elements corresponding to
             all grid points inside the cloudbox plus the cloud_box_limits*/

      p2gridpos(p_gp,
                p_grid[Range(cloudbox_limits[0],
                             2,
                             (cloudbox_limits[1] - cloudbox_limits[0]))],
                p_grid[Range(cloudbox_limits[0],
                             (cloudbox_limits[1] - cloudbox_limits[0]) + 1)]);
      gridpos(lat_gp,
              lat_grid[Range(cloudbox_limits[2],
                             2,
                             (cloudbox_limits[3] - cloudbox_limits[2]))],
              lat_grid[Range(cloudbox_limits[2],
                             (cloudbox_limits[3] - cloudbox_limits[2]) + 1)]);
      gridpos(lon_gp,
              lon_grid[Range(cloudbox_limits[4],
                             2,
                             (cloudbox_limits[5] - cloudbox_limits[4]))],
              lon_grid[Range(cloudbox_limits[4],
                             (cloudbox_limits[5] - cloudbox_limits[4]) + 1)]);

      //interpolation weights corresponding to pressure, latitude and
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1] - cloudbox_limits[0]) + 1, 2);
      Matrix itw_lat((cloudbox_limits[3] - cloudbox_limits[2]) + 1, 2);
      Matrix itw_lon((cloudbox_limits[5] - cloudbox_limits[4]) + 1, 2);

      interpweights(itw_p, p_gp);
      interpweights(itw_lat, lat_gp);
      interpweights(itw_lon, lon_gp);

      // interpolation - pressure grid
      for (Index lat_index = 0;
           lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
           ++lat_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = doit_i_field(f_index,
                                                       Range(joker),
                                                       lat_index,
                                                       lon_index,
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_p(
                    Range(joker), lat_index, lon_index, za_index, aa_index, i);

                interp(target_field, itw_p, source_field, p_gp);
              }
            }
          }
        }
      }
      //interpolation latitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = doit_i_field(f_index,
                                                       p_index,
                                                       Range(joker),
                                                       lon_index,
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_lat(
                    p_index, Range(joker), lon_index, za_index, aa_index, i);

                interp(target_field, itw_lat, source_field, lat_gp);
              }
            }
          }
        }
      }
      //interpolation -longitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lat_index = 0;
             lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
             ++lat_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = doit_i_field(f_index,
                                                       p_index,
                                                       lat_index,
                                                       Range(joker),
                                                       za_index,
                                                       aa_index,
                                                       i);

                ConstVectorView source_field = scat_i_lon(
                    p_index, lat_index, Range(joker), za_index, aa_index, i);

                interp(target_field, itw_lon, source_field, lon_gp);
              }
            }
          }
        }
      }  //end of interpolation
    }    // end of frequency loop
  }      //ends atmosphere_dim = 3
}



